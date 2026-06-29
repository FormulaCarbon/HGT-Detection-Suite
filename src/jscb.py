import argparse
from pathlib import Path

from jscbhelpers import *

def main(genbank: Path, output: Path, verbose: bool, gpu: bool, gpu_batch_size: int) -> int:
	output.mkdir(exist_ok=True, parents=True)
	script_path =  Path(__file__).parent.resolve()
 
	rRNA_iterator = SeqIO.parse(script_path / 'E_coli_16S.gb', 'genbank')
	rRNA_record = next(rRNA_iterator)
	combined_record = None
	CDS_locus_tags = list() # Ordered list of CDS locus_tags so that we can assign them back to gene numbers later
	for seq_record in SeqIO.parse(genbank, 'genbank'):
		for feature in seq_record.features:
			if feature.type == 'CDS':
				locus_tag = feature.qualifiers['locus_tag'][0]
				CDS_locus_tags.append(locus_tag)
		if combined_record is None:
			combined_record = seq_record
		else:
			combined_record = combined_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + seq_record
	combined_gbk_path = output / 'combined.gbk'
	if combined_record is None:
		print("combined_record is None")
		return 1
	SeqIO.write(combined_record, combined_gbk_path, 'genbank')
	genome, genes_raw, rrna = load_from_genbank(combined_gbk_path)
	genes_all = genes_raw

	# The Fortran filters out short genes by:
	# lena = end - start - 2; if (lena-3 > mingenelen) keep, with mingenelen=3
	# => keep if lena > 6, i.e. (end-start-2) > 6 => length_nt > 8
	# We'll implement the same.
	mingenelen = 3
	genes_used = []
	for g in genes_all:
		lena = g.end - g.start - 2
		if (lena - 3) > mingenelen:
			genes_used.append(g)

	if not genes_used:
		print("No genes passed filters (CDS length % 3 == 0 and minimum length).")
		return 1

	iy = genome_to_iy(genome)
	iaa, iab = build_iaa_iab()

	ncodonfreq, naminofreq, gene_len_internal = compute_gene_codonfreq_and_aminofreq(
		genes=genes_used, iy=iy, iaa=iaa, iab=iab
	)

	gene_to_cluster_used, clusters, _, _, _ = cluster_jscb(
		ncodonfreq=ncodonfreq,
		naminofreq=naminofreq,
		gene_len_internal=gene_len_internal,
		iaa=iaa,
		iab=iab,
		thres=0.8,
		thres2=0.0,
		thres3=0.995,
		clusthres=1.0,
		use_gpu=gpu,
		gpu_batch_size=max(1, gpu_batch_size),
		verbose=verbose,
	)

	# map clustered genes_used back to genes_all by exact coordinates+strand
	cluster_by_coord = {}
	for idx_used, g in enumerate(genes_used, start=1):
		cluster_by_coord[(g.start, g.end, g.strand)] = gene_to_cluster_used[idx_used]

	gene_to_cluster_all = [0] * (len(genes_all) + 1)
	for idx_all, g in enumerate(genes_all, start=1):
		gene_to_cluster_all[idx_all] = cluster_by_coord.get((g.start, g.end, g.strand), 0)

	cluster_sizes = {}
	for idx_used in range(1, len(gene_to_cluster_used)):
		cid = gene_to_cluster_used[idx_used]
		if cid > 0:
			cluster_sizes[cid] = cluster_sizes.get(cid, 0) + 1

	# match Fortran lena(m) = end - start - 2 in genes_all index space
	gene_lena = [0] + [(g.end - g.start - 2) for g in genes_all]

	if verbose:
		mapped_nonzero = sum(1 for cid in gene_to_cluster_all[1:] if cid > 0)
		unmapped = len(genes_all) - mapped_nonzero
		native_cid = 0
		native_size = 0
		for cid, csize in cluster_sizes.items():
			if csize > native_size:
				native_cid = cid
				native_size = csize
		print(
			f"DEBUG stats: genes_all={len(genes_all)} genes_used={len(genes_used)} "
			f"mapped={mapped_nonzero} unmapped={unmapped} nclusters={len(clusters)} "
			f"native_cluster={native_cid} native_size={native_size}",
		)

	write_output_gi(output / "JSCB_output.gi", gene_to_cluster_all, cluster_sizes, gene_lena)
	call_genomic_islands(output / "JSCB_output.tsv", gene_to_cluster_all, genes_all, rrna)
	
	# Post Processing
	jscb_output_path = output / 'JSCB_output.tsv'
	genomic_islands = load_genomic_island_ranges(jscb_output_path)
	genomic_island_genes = build_genomic_island_gene_map(combined_gbk_path, genomic_islands)
	genomic_island_coordinates, genomic_island_genes_by_contig = build_genomic_island_contig_details(
		genbank,
		genomic_island_genes,
	)

	# Now we make a human-readable output table
	output_table_path = output / 'genomic_islands_summary.tsv'
	output_table = open(output_table_path, 'w')
	output_table.write('GI_ID\tcontig\tstart\tend\tgenes\n')

	for gi_name in genomic_island_coordinates:
		if len(genomic_island_coordinates[gi_name]) > 1:
			# This means the GI spans multiple contigs
			# Because the order of contigs in the original file is likely random,
			# This means in reality there are two GIs
			GI_counter = 0
			for contig_name in genomic_island_coordinates[gi_name]:
				GI_counter += 1
				output_gi_name = gi_name + '_' + str(GI_counter)
				gene_string = ','.join(genomic_island_genes_by_contig[gi_name][contig_name])
				output_string = '\t'.join([ output_gi_name, contig_name, str(genomic_island_coordinates[gi_name][contig_name][0]), str(genomic_island_coordinates[gi_name][contig_name][1]), gene_string ])
				output_table.write(output_string + '\n')
		else:
			for contig_name in genomic_island_coordinates[gi_name]:
				gene_string = ','.join(genomic_island_genes_by_contig[gi_name][contig_name])
				output_string = '\t'.join([ gi_name, contig_name, str(genomic_island_coordinates[gi_name][contig_name][0]), str(genomic_island_coordinates[gi_name][contig_name][1]), gene_string ])
				output_table.write(output_string + '\n')
	output_table.close()

	# We also make a human readable version of JSCB_output.gi
	output_table_path = output / 'hgt_genes_summary.tsv'
	output_table = open(output_table_path, 'w')
	output_table.write('locus_tag\tcluster_id\tcluster_size\tgene_length\n')
	jscb_output_path = output / 'JSCB_output.gi'

	with open(jscb_output_path) as jscb_output:
		jscb_output.readline() # throw away header
		for line in jscb_output:
			line_list = line.rstrip().split()
			gene_index = int(line_list[0]) - 1

			locus_tag = CDS_locus_tags[gene_index]

			output_string = '\t'.join([locus_tag] + line_list[1:])
			output_table.write(output_string + '\n')

	output_table.close()


	return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank", type=Path, help="Path to annotated genbank")
    parser.add_argument("-o", "--output", type=Path, default="/", help="Path to output directory")
    parser.add_argument("-v", "--verbose", action='store_true', help="Enable verbose logging")
    parser.add_argument("--gpu", action='store_true', help="Enable GPU-accelerated merge-step divergence checks when CuPy/CUDA is available")
    parser.add_argument("--gpu-batch-size", type=int, default=64, help="Number of cluster-pairs to process per GPU batch during merge step")
    
    args = parser.parse_args()
    
    main(args.genbank, args.output, args.verbose, args.gpu, args.gpu_batch_size)