import math
import random
from pathlib import Path

import numpy as np
from scipy.special import gammainc
from Bio import SeqIO

try:
	import cupy as cp
except Exception:
	cp = None
 
from dataclasses import dataclass

@dataclass
class Gene:
	strand: str  # '+' or '-'
	start: int   # 1-based inclusive
	end: int	 # 1-based inclusive

	@property
	def length_nt(self) -> int:
		return self.end - self.start + 1

BASE_TO_IY = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
COMP_IY = {1: 2, 2: 1, 3: 4, 4: 3}
INV_LOG2 = 1.0 / math.log(2.0)

def get_cupy_module():
	if cp is None:
		return None
	try:
		if cp.cuda.runtime.getDeviceCount() < 1:
			return None
	except Exception:
		return None
	return cp

def build_synonymous_codon_groups(iaa: list[list[int]], iab: list[int]) -> list[np.ndarray]:
	return [np.array( iaa[aa_idx][:iab[aa_idx]], dtype=np.int32 ) for aa_idx in range(20)]

def sig(rnn_array, djsmx_array):
	snn = rnn_array * math.log(2.0) * djsmx_array
	return gammainc(30.0, snn)

def build_iaa_iab() -> tuple:
	"""
	Returns:
	  iaa: 20 x 6 (1-based codonIndex values), 0 where unused
	  iab: length 20, count of synonymous codons in that amino acid group
	"""
	iaa = [[0] * 6 for _ in range(20)]
	def set_row(r: int, vals: list[int]) -> None:
		for c, v in enumerate(vals):
			iaa[r - 1][c] = v

	set_row(1,  [42, 43])
	set_row(2,  [41, 44, 58, 59, 57, 60])
	set_row(3,  [26, 27, 25])
	set_row(4,  [28])
	set_row(5,  [74, 75, 73, 76])
	set_row(6,  [46, 47, 45, 48, 34, 35])
	set_row(7,  [62, 63, 61, 64])
	set_row(8,  [30, 31, 29, 32])
	set_row(9,  [78, 79, 77, 80])
	set_row(10, [38, 39])
	set_row(11, [54, 55])
	set_row(12, [53, 56])
	set_row(13, [22, 23])
	set_row(14, [21, 24])
	set_row(15, [70, 71])
	set_row(16, [69, 72])
	set_row(17, [50, 51])
	set_row(18, [52])
	set_row(19, [66, 67, 65, 68, 33, 36])
	set_row(20, [82, 83, 81, 84])

	iab = [
		2, 6, 3, 1, 4, 6, 4, 4, 4, 2,
		2, 2, 2, 2, 2, 2, 2, 1, 6, 4
	]
	return iaa, iab

def genome_to_iy(genome: str):
	iy = np.zeros(len(genome) + 1, dtype=np.int8)
	for i, ch in enumerate(genome, start=1):
		try:
			iy[i] = BASE_TO_IY[ch]
		except KeyError:
			raise ValueError(f"problem with genome file at position {i}: {ch!r}")
	return iy

def codon_index_from_iy(k1: int, k2: int, k3: int) -> int:
	return (16 * k1) + (4 * k2) + k3

def make_codon_to_aa(iaa, iab, ih):
	codon_to_aa = np.zeros(ih+1, dtype=np.int16)

	for aa in range(20):
		for cod in iaa[aa][:iab[aa]]:
			codon_to_aa[cod] = aa + 1  # 1-based AA group
	return codon_to_aa

def compute_gene_codonfreq_and_aminofreq(
	genes: list[Gene],
	iy,
	iaa: list[list[int]],
	iab: list[int],
	ih0: int = 21,
	ih: int = 84,
):
    
	codon_to_aa = make_codon_to_aa(iaa, iab, ih)
	aa_codons = [np.where(codon_to_aa == aa + 1)[0] for aa in range(20)]
	"""
	Returns:
	  ncodonfreq[geneIndex][codonIndex] (1-based gene index; codonIndex 0..ih)
	  naminofreq[geneIndex][codonIndex]
	  gene_len_internal[geneIndex] = irend - lend - 2 (Fortran's len(j))
	"""


	ngenea = len(genes)
	# allocate with padding for 1-based gene indexing and codon index up to ih
	ncodonfreq = np.zeros((ngenea + 1, ih + 1), dtype=np.int32)
	naminofreq = np.zeros((ngenea + 1, ih + 1), dtype=np.int32)
	gene_len_internal = np.zeros(ngenea + 1, dtype=np.int32)

	# init
	ncodonfreq[:, ih0:ih+1] = 1

	for g, gene in enumerate(genes, start=1):
		lend = gene.start
		irend = gene.end
		gene_len_internal[g] = irend - lend - 2

		ncodons = gene_len_internal[g] // 3

		if ncodons <= 0:
			continue
  
		if gene.strand == '+':
			for ii in range(ncodons):
				i = lend + ii * 3

				c1 = iy[i]
				c2 = iy[i + 1]
				c3 = iy[i + 2]

				codon = (c1 << 4) + (c2 << 2) + c3
				ncodonfreq[g, codon] += 1

		else:
			for ii in range(ncodons):
				i = irend - ii * 3

				c1 = COMP_IY[iy[i]]
				c2 = COMP_IY[iy[i - 1]]
				c3 = COMP_IY[iy[i - 2]]

				codon = (c1 << 4) + (c2 << 2) + c3
				ncodonfreq[g, codon] += 1

		# Build amino frequency once codon counts are complete
		for codons in aa_codons:
			if codons.size > 0:
				naminofreq[g, codons] = int(np.sum(ncodonfreq[g, codons]))

	return ncodonfreq, naminofreq, gene_len_internal

def load_genomic_island_ranges(path: Path) -> dict[str, tuple[int, int]]:
	genomic_islands = {}
	with open(path, encoding="utf-8") as jscb_output:
		for i, line in enumerate(jscb_output):
			if i == 0:
				continue
			line_list = line.rstrip().split('\t')
			if len(line_list) < 3:
				continue
			genomic_islands[line_list[0]] = (int(line_list[1]), int(line_list[2]))
	return genomic_islands

def build_genomic_island_gene_map(
	combined_gbk_path: Path,
	genomic_islands: dict[str, tuple[int, int]],
) -> dict[str, set[str]]:
    
	genomic_island_genes = {gi_name: set() for gi_name in genomic_islands}
	for seq_record in SeqIO.parse(combined_gbk_path, 'genbank'):
		for feature in seq_record.features:
			if feature.type != 'gene':
				continue
			locus_tags = feature.qualifiers.get('locus_tag')
			if not locus_tags:
				continue
			gene_start = int(feature.location.start)
			gene_end = int(feature.location.end)
			locus_tag = locus_tags[0]
			for gi_name, (gi_start, gi_end) in genomic_islands.items():
				if gene_start <= gi_end and gene_end >= gi_start:
					genomic_island_genes[gi_name].add(locus_tag)
	return genomic_island_genes

def build_genomic_island_contig_details(
	genbank_path: Path,
	genomic_island_genes: dict[str, set[str]],
) -> tuple:
    
	genomic_island_coordinates = {gi_name: {} for gi_name in genomic_island_genes}
	genomic_island_genes_by_contig = {gi_name: {} for gi_name in genomic_island_genes}

	locus_to_islands = {}
	for gi_name, locus_tags in genomic_island_genes.items():
		for locus_tag in locus_tags:
			locus_to_islands.setdefault(locus_tag, []).append(gi_name)

	for seq_record in SeqIO.parse(genbank_path, 'genbank'):
		contig_name = seq_record.id
		for feature in seq_record.features:
			if feature.type != 'gene':
				continue
			locus_tags = feature.qualifiers.get('locus_tag')
			if not locus_tags:
				continue
			locus_tag = locus_tags[0]
			gi_names = locus_to_islands.get(locus_tag)
			if not gi_names:
				continue
			start_coordinate = int(feature.location.start)
			end_coordinate = int(feature.location.end)
			for gi_name in gi_names:
				if contig_name in genomic_island_coordinates[gi_name]:
					if start_coordinate < genomic_island_coordinates[gi_name][contig_name][0]:
						genomic_island_coordinates[gi_name][contig_name][0] = start_coordinate
					if end_coordinate > genomic_island_coordinates[gi_name][contig_name][1]:
						genomic_island_coordinates[gi_name][contig_name][1] = end_coordinate
					genomic_island_genes_by_contig[gi_name][contig_name].append(locus_tag)
				else:
					genomic_island_coordinates[gi_name][contig_name] = [start_coordinate, end_coordinate]
					genomic_island_genes_by_contig[gi_name][contig_name] = [locus_tag]

	return genomic_island_coordinates, genomic_island_genes_by_contig

def load_from_genbank(path: Path):
	genes = []
	rrna = []

	recs = list(SeqIO.parse(path, "genbank"))
	if not recs:
		raise ValueError("No GenBank records found")

	rec = recs[0]
	genome = sanitize_genome(str(rec.seq))

	for feature in rec.features:
		start = int(feature.location.start) + 1
		end = int(feature.location.end)
		s = feature.location.strand
		if s == 1:
			strand = "+"
		elif s == -1:
			strand = "-"
		else:
			strand = "+"
		if feature.type == "CDS":
			if "join" not in str(feature.location):
				genes.append(Gene(strand=strand, start=start, end=end))
		elif feature.type == "rRNA":
			rrna.append((start, end))
   
	filtered_genes = [g for g in genes if (g.length_nt % 3) == 0]
	return genome, filtered_genes, rrna

def sanitize_genome(seq: str) -> str:
	nucs = ['A', 'T', 'C', 'G']
	out = []
	for ch in seq.upper():
		if ch in BASE_TO_IY:
			out.append(ch)
		else:
			out.append(random.choice(nucs))
	return "".join(out)

def _divgn_cluster_vs_many_gpu(
	xp,
	c1_codon: np.ndarray,
	c1_amino: np.ndarray,
	c1_len: int,
	c2_codon_batch: np.ndarray,
	c2_amino_batch: np.ndarray,
	c2_len_batch: np.ndarray,
	synonymous_codons: list[np.ndarray],
) -> np.ndarray:
	c1_codon_xp = xp.asarray(c1_codon, dtype=xp.float64)
	c1_amino_xp = xp.asarray(c1_amino, dtype=xp.float64)
	c2_codon_xp = xp.asarray(c2_codon_batch, dtype=xp.float64)
	c2_amino_xp = xp.asarray(c2_amino_batch, dtype=xp.float64)
	c2_len_xp = xp.asarray(c2_len_batch, dtype=xp.float64)
	c1_len_f = float(c1_len)
	lengtot_xp = c1_len_f + c2_len_xp

	ent1 = xp.float64(0.0)
	ent2 = xp.zeros(c2_len_xp.shape[0], dtype=xp.float64)
	enttot = xp.zeros(c2_len_xp.shape[0], dtype=xp.float64)

	for codons_np in synonymous_codons:
		codons = xp.asarray(codons_np)
		first = int(codons_np[0])

		c1_den = c1_amino_xp[codons]
		c1_pp = c1_codon_xp[codons] / c1_den
		ent1 += xp.sum(c1_pp * xp.log(c1_pp)) * (c1_amino_xp[first] / c1_len_f)

		c2_den = c2_amino_xp[:, codons]
		c2_pp = c2_codon_xp[:, codons] / c2_den
		ent2 += xp.sum(c2_pp * xp.log(c2_pp), axis=1) * (c2_amino_xp[:, first] / c2_len_xp)

		denom_tot = c1_amino_xp[first] + c2_amino_xp[:, first]
		pp_tot = (c1_codon_xp[codons][None, :] + c2_codon_xp[:, codons]) / denom_tot[:, None]
		enttot += xp.sum(pp_tot * xp.log(pp_tot), axis=1) * (denom_tot / lengtot_xp)

	divergence_xp = INV_LOG2 * ((-enttot) - ((c1_len_f * (-ent1) + c2_len_xp * (-ent2)) / lengtot_xp))
	if xp is np:
		return divergence_xp
	return xp.asnumpy(divergence_xp)

def cluster_jscb(
    ncodonfreq: np.ndarray,
    naminofreq: np.ndarray,
    gene_len_internal,
    iaa: list[list[int]],
    iab: list[int],
    ih0: int = 21,
    ih: int = 84,
    thres: float = 0.8,
    thres2: float = 0.0,
    thres3: float = 0.995,
    clusthres: float = 1.0,
	use_gpu: bool = False,
	gpu_batch_size: int = 64,
	verbose: bool = False,
):

    ngenes = len(gene_len_internal) - 1
    max_clusters = ngenes + 1

    # -------------------------
    # Precompute AA mapping
    # -------------------------
    aa_first = np.array([iaa[i][0] for i in range(20)], dtype=np.int32)
    synonymous_codons = build_synonymous_codon_groups(iaa, iab)
    gpu_xp = get_cupy_module() if use_gpu else None
    gpu_enabled = gpu_xp is not None
    if use_gpu and (not gpu_enabled) and verbose:
        print("GPU requested but unavailable; falling back to CPU.")

    # -------------------------
    # Cluster state (FAST)
    # -------------------------
    cluster_codon = np.zeros((max_clusters, ih + 1), dtype=np.int32)
    cluster_amino = np.zeros((max_clusters, ih + 1), dtype=np.int32)
    cluster_len = np.zeros(max_clusters, dtype=np.int32)

    gene_to_cluster = np.zeros(ngenes + 1, dtype=np.int32)
    clusters = []

    # -------------------------
    # Init first cluster
    # -------------------------
    def init_cluster(cid: int, gene_idx: int):
        clusters.append([gene_idx])
        gene_to_cluster[gene_idx] = cid

        cluster_codon[cid, ih0:ih + 1] = ncodonfreq[gene_idx, ih0:ih + 1]
        cluster_amino[cid, ih0:ih + 1] = naminofreq[gene_idx, ih0:ih + 1]

        cluster_len[cid] = cluster_amino[cid, aa_first].sum()

    def add_gene(cid: int, gene_idx: int):
        clusters[cid - 1].append(gene_idx)
        gene_to_cluster[gene_idx] = cid

        cluster_codon[cid, ih0:ih + 1] += ncodonfreq[gene_idx, ih0:ih + 1]
        cluster_amino[cid, ih0:ih + 1] += naminofreq[gene_idx, ih0:ih + 1]

        cluster_len[cid] += gene_len_internal[gene_idx] // 3

    # -------------------------
    # Sequential clustering
    # -------------------------
    init_cluster(1, 1)
    current_cid = 1

    for g in range(2, ngenes + 1):
        gene_len_codons = gene_len_internal[g] // 3
        lengtot = cluster_len[current_cid] + gene_len_codons

        ddiv = divgn2_gene_cluster(
            ncodonfreq[g], naminofreq[g], gene_len_codons,
            cluster_codon[current_cid], cluster_amino[current_cid],
            cluster_len[current_cid],
            lengtot, iaa, iab
        )

        pbc = sig(float(lengtot), abs(ddiv))

        if pbc < thres:
            add_gene(current_cid, g)
        else:
            current_cid += 1
            init_cluster(current_cid, g)

    nclus = len(clusters)

    # -------------------------
    # Merge step (OPTIMIZED)
    # -------------------------
    cluster_entropy = np.full(max_clusters, np.nan, dtype=np.float64)

    changed = True
    while changed:
        changed = False

        for i in range(1, nclus):
            if cluster_len[i] == 0:
                continue
            candidate_js = [j for j in range(i + 1, nclus + 1) if cluster_len[j] > 0]
            if not candidate_js:
                continue

            merged_with_j = None

            if gpu_enabled and len(candidate_js) >= 2:
                for start in range(0, len(candidate_js), max(1, gpu_batch_size)):
                    batch_js = candidate_js[start:start + max(1, gpu_batch_size)]
                    djs_vals = _divgn_cluster_vs_many_gpu(
                        gpu_xp,
                        cluster_codon[i],
                        cluster_amino[i],
                        int(cluster_len[i]),
                        cluster_codon[batch_js],
                        cluster_amino[batch_js],
                        cluster_len[batch_js],
                        synonymous_codons,
                    )
                    for idx, j in enumerate(batch_js):
                        lengtot = cluster_len[i] + cluster_len[j]
                        pbc = sig(float(lengtot), float(djs_vals[idx]))
                        if pbc < thres3:
                            merged_with_j = j
                            break
                    if merged_with_j is not None:
                        break
            else:
                for j in candidate_js:
                    lengtot = cluster_len[i] + cluster_len[j]

                    ent_i = cluster_entropy[i]
                    if np.isnan(ent_i):
                        ent_i = rentro_cluster(cluster_codon[i], cluster_amino[i], cluster_len[i], iaa, iab)
                        cluster_entropy[i] = ent_i

                    ent_j = cluster_entropy[j]
                    if np.isnan(ent_j):
                        ent_j = rentro_cluster(cluster_codon[j], cluster_amino[j], cluster_len[j], iaa, iab)
                        cluster_entropy[j] = ent_j

                    djs = INV_LOG2 * (
                        rentrotot_cluster_cluster(
                            cluster_codon[i], cluster_amino[i],
                            cluster_codon[j], cluster_amino[j],
                            lengtot, iaa, iab,
                        )
                        - ((cluster_len[i] * ent_i + cluster_len[j] * ent_j) / float(lengtot))
                    )
                    pbc = sig(float(lengtot), djs)
                    if pbc < thres3:
                        merged_with_j = j
                        break

            if merged_with_j is not None:
                j = merged_with_j

                # merge j -> i (FAST vectorized update)
                clusters[i - 1].extend(clusters[j - 1])

                cluster_codon[i] += cluster_codon[j]
                cluster_amino[i] += cluster_amino[j]
                cluster_len[i] += cluster_len[j]

                # invalidate j (NO deletion)
                cluster_codon[j].fill(0)
                cluster_amino[j].fill(0)
                cluster_len[j] = 0
                clusters[j - 1] = []
                cluster_entropy[i] = np.nan
                cluster_entropy[j] = np.nan

                for g in clusters[i - 1]:
                    gene_to_cluster[g] = i

                changed = True
                break

            if changed:
                break

    # compress cluster IDs (remove empty ones)
    mapping = {}
    new_clusters = []
    new_codon = []
    new_amino = []
    new_len = []

    cid_new = 0
    for cid in range(1, nclus + 1):
        if len(clusters[cid - 1]) == 0:
            continue

        cid_new += 1
        mapping[cid] = cid_new

        new_clusters.append(clusters[cid - 1])
        new_codon.append(cluster_codon[cid])
        new_amino.append(cluster_amino[cid])
        new_len.append(cluster_len[cid])

    clusters = new_clusters
    cluster_codon = np.array(new_codon)
    cluster_amino = np.array(new_amino)
    cluster_len = np.array(new_len)

    # remap gene_to_cluster
    for g in range(1, ngenes + 1):
        if gene_to_cluster[g] in mapping:
            gene_to_cluster[g] = mapping[gene_to_cluster[g]]

    while True:
        nclus = len(clusters)

        rpos = np.zeros((nclus + 1, nclus + 1), dtype=np.float32)

        for cid, genes in enumerate(clusters, start=1):
            for g in genes:
                if g <= 1 or g >= ngenes:
                    continue

                c_left = gene_to_cluster[g - 1]
                c_right = gene_to_cluster[g + 1]

                if c_left == c_right and c_left != 0:
                    rpos[cid, c_left] += 1

        merged = False

        for i in range(1, nclus):
            for j in range(i + 1, nclus + 1):

                rp = (rpos[i, j] / max(1, len(clusters[i - 1])) +
                      rpos[j, i] / max(1, len(clusters[j - 1]))) / 2.0

                if rp > clusthres:

                    clusters[i - 1].extend(clusters[j - 1])

                    cluster_codon[i] += cluster_codon[j]
                    cluster_amino[i] += cluster_amino[j]
                    cluster_len[i] += cluster_len[j]

                    clusters[j - 1] = []

                    for g in clusters[i - 1]:
                        gene_to_cluster[g] = i

                    merged = True
                    break

            if merged:
                break

        if not merged:
            break

    return gene_to_cluster, clusters, cluster_codon, cluster_amino, cluster_len

def divgn2_gene_cluster(
	gene_codon: list[int], gene_am: list[int], gene_len_codons: int,
	cl_codon: list[int], cl_am: list[int], icluslen: int,
	lengtot: int,
	iaa: list[list[int]], iab: list[int],
) -> float:

	hhh = 2.0
	xxx = 1.0 / math.log(hhh)
	return xxx * (
		rentrotot_gene_cluster(gene_codon, gene_am, cl_codon, cl_am, lengtot, iaa, iab)
		- ((gene_len_codons) * rentro_gene(gene_codon, gene_am, gene_len_codons, iaa, iab)
		   + icluslen * rentro_cluster(cl_codon, cl_am, icluslen, iaa, iab)) / float(lengtot)
	)

def rentrotot_gene_cluster(
	gene_codon: list[int],
	gene_am: list[int],
	cl_codon: list[int],
	cl_am: list[int],
	lengtot: int,
	iaa: list[list[int]],
	iab: list[int],
) -> float:
	# Fortran rentrotot2(id1,id2,lengtot) (id1 is gene, id2 is cluster)
	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		denom = float(gene_am[iaa[aa_idx][0]] + cl_am[iaa[aa_idx][0]])
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = (gene_codon[cod] + cl_codon[cod]) / denom
			rentroi += pp * math.log(pp)
		pp0 = (gene_am[iaa[aa_idx][0]] + cl_am[iaa[aa_idx][0]]) / float(lengtot)
		rentroi *= pp0
		val += rentroi
	return -val


def divgn_cluster_cluster(
	c1_codon: list[int], c1_am: list[int], icluslen1: int,
	c2_codon: list[int], c2_am: list[int], icluslen2: int,
	lengtot: int,
	iaa: list[list[int]], iab: list[int],
) -> float:
	
	hhh = 2.0
	xxx = 1.0 / math.log(hhh)
	return xxx * (
		rentrotot_cluster_cluster(c1_codon, c1_am, c2_codon, c2_am, lengtot, iaa, iab)
		- (icluslen1 * rentro_cluster(c1_codon, c1_am, icluslen1, iaa, iab)
		   + icluslen2 * rentro_cluster(c2_codon, c2_am, icluslen2, iaa, iab)) / float(lengtot)
	)
 
def rentro_cluster(
	cluster_codon: list[int],
	cluster_am: list[int],
	icluslen: int,
	iaa: list[list[int]],
	iab: list[int],
) -> float:

	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = cluster_codon[cod] / float(cluster_am[cod])
			rentroi += pp * math.log(pp)
		rentroi *= (cluster_am[iaa[aa_idx][0]] / float(icluslen))
		val += rentroi
	return -val


def rentro_gene(
	gene_codon: list[int],
	gene_am: list[int],
	gene_len_codons: int,
	iaa: list[list[int]],
	iab: list[int],
) -> float:

	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = gene_codon[cod] / float(gene_am[cod])
			rentroi += pp * math.log(pp)
		rentroi *= (gene_am[iaa[aa_idx][0]] / float(gene_len_codons))
		val += rentroi
	return -val


def rentrotot_cluster_cluster(
	c1_codon: list[int],
	c1_am: list[int],
	c2_codon: list[int],
	c2_am: list[int],
	lengtot: int,
	iaa: list[list[int]],
	iab: list[int],
) -> float:

	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		denom = float(c1_am[iaa[aa_idx][0]] + c2_am[iaa[aa_idx][0]])
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = (c1_codon[cod] + c2_codon[cod]) / denom
			rentroi += pp * math.log(pp)
		pp0 = (c1_am[iaa[aa_idx][0]] + c2_am[iaa[aa_idx][0]]) / float(lengtot)
		rentroi *= pp0
		val += rentroi
	return -val

def write_output_gi(
	out_path: Path,
	gene_to_cluster: list[int],
	cluster_sizes: dict[int, int],
	gene_lena: list[int],
):
	# geneIndex clusterId clusterSize geneLength
	with open(out_path, "w", encoding="utf-8") as f:
		f.write(f"gene_number\tclus_id\tclus_size\tgene_length\n")
		for g in range(1, len(gene_to_cluster)):
			cid = gene_to_cluster[g]
			csize = cluster_sizes.get(cid, 0) if cid > 0 else 0
			f.write(f"{g}\t{cid}\t{csize}\t{gene_lena[g]}\n")


def call_genomic_islands(
	out_tsv: Path,
	gene_to_cluster: list[int],
	genes: list[Gene],
	rrna: list[tuple[int, int]],
	min_genes: int = 8,
):
	# Find native cluster = max cluster size (by gene count)
	counts = {}
	for g in range(1, len(gene_to_cluster)):
		cid = gene_to_cluster[g]
		if cid > 0:
			counts[cid] = counts.get(cid, 0) + 1
	if not counts:
		with open(out_tsv, "w", encoding="utf-8") as outfile:
			outfile.write("No Genomic Islands Identified\n")
		return
	native_cluster = max(counts, key=counts.get)

	actstart = [0] + [g.start for g in genes]  # 1-based
	actend = [0] + [g.end for g in genes]	  # 1-based

	prv = -1
	gmode = 0
	gicounter = 0

	with open(out_tsv, "w", encoding="utf-8") as outfile:
		for gene_idx in range(1, len(gene_to_cluster)):
			cid = gene_to_cluster[gene_idx]
			if cid != native_cluster:
				if gene_idx != prv + 1:
					startgeneno = gene_idx
					startcoord = actstart[gene_idx]
					gmode = 1
				prv = gene_idx
			elif gmode == 1:
				endgeneno = gene_idx
				gilength = endgeneno - startgeneno
				endcoord = actend[gene_idx - 1]

				if gilength >= min_genes:
					rnais = 0
					for x, y in rrna:
						# exact same containment test as jscb.py
						if startcoord < x and startcoord < y and endcoord > x and endcoord > y:
							rnais = 1
							break
					if rnais == 0:
						gicounter += 1
						if gicounter == 1:
							outfile.write("clus_id\tstart\tend\n")
						outfile.write(f"GI-{gicounter}\t{startcoord}\t{actend[gene_idx - 1]}\n")

				gmode = 0

		if gicounter == 0:
			outfile.write("No Genomic Islands Identified\n")