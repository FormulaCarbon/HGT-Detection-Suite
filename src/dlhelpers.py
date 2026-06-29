from pathlib import Path
import subprocess
import urllib.request
import time
import shutil
import gzip

from iohelpers import header, error
                
def download_dataset(accessions_file: Path, dehydrated_file: Path, db_path: Path, max_workers = 10, api_key = None):
    header("Downloading Dehydrated Dataset...")
    subprocess.run([
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        accessions_file,
        "--dehydrated",
        "--filename",
        dehydrated_file,  
        "--include",
        "protein",
        ] + (["--api-key", api_key] if api_key is not None else []),
        shell = True
    )
    print()
    
    header("Unzipping Dehydrated Dataset...")
    subprocess.run([
        "tar",
        "-xf",
        dehydrated_file,
        "-C",
        db_path,
        ], 
        shell = True
    )
    print()
    
    header("Rehydrating Dataset...")
    subprocess.run([
        "datasets",
        "rehydrate",
        "--directory",
        db_path,
        "--max-workers",
        str(max_workers)
        ] + (["--api-key", api_key] if api_key is not None else [])
    )
    print()
    
def download_url(url: str, out_path: Path, retries: int = 3, wait: float = 1.0) -> int:

    out_path.parent.mkdir(parents = True, exist_ok = True)
    for i in range(retries):
        try:
            with urllib.request.urlopen(url) as response:
                data = response.read()

                with open(out_path, 'wb') as file:
                    file.write(data)
            return 0
        except:
            time.sleep(wait)
    
    return 1

def unzip_fastas(targetPath: Path):
    for gz in targetPath.glob("*.faa.gz"): 
        out = Path(str(gz)[:-3])  # remove .gz 
        if out.exists(): 
            continue 
        print("unzipping", gz.name) 
        with gzip.open(gz, "rb") as f_in, open(out, "wb") as f_out: 
            shutil.copyfileobj(f_in, f_out) 
    