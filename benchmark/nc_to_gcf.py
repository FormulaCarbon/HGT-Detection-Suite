import subprocess
import json


def get_gcf_accession(nc_accession):
    search = subprocess.Popen(
        ['esearch', '-db', 'assembly', '-query', nc_accession],
        stdout = subprocess.PIPE
    )
    summary = subprocess.Popen(
        ['esummary'],
        stdin = search.stdout,
        stdout = subprocess.PIPE
    )
    xtract = subprocess.Popen(
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'RefSeq'],
        stdin = summary.stdout,
        stdout = subprocess.PIPE,
        text=True
    )
    
    output, _ = xtract.communicate() 
    return output

# Example usage
if __name__ == "__main__":
    ids = {}
    with open("benchmark/nc_accessions.txt", 'r') as f:
        for line in f:
            nc_id = line.strip()
            gcf_id = get_gcf_accession(nc_id).strip()
            ids[nc_id] = gcf_id
            print(nc_id, gcf_id)
            
    with open("benchmark/accessions.txt", 'w') as out:
        json.dump(ids, out)