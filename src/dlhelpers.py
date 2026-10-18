from pathlib import Path
import subprocess
import urllib.request
import time
import shutil
import gzip
import os

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
        "genome,protein,gbff",
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

def split_by_filetype(
    source_dir: Path,
    faa_dir: Path,
    fna_dir: Path,
    gbff_dir: Path,
    copy: bool = False,
) -> dict[str, int]:
    """Sort .faa, .fna, and .gbff files from source_dir into separate folders."""
    destinations = {
        ".faa": Path(faa_dir),
        ".fna": Path(fna_dir),
        ".gbff": Path(gbff_dir),
    }

    for dest in destinations.values():
        dest.mkdir(parents=True, exist_ok=True)

    transfer = shutil.copy2 if copy else shutil.move
    counts = {ext: 0 for ext in destinations}

    # scandir iterates lazily, which is better than building a list of ~180k paths
    with os.scandir(source_dir) as entries:
        for entry in entries:
            if not entry.is_file():
                continue
            ext = Path(entry.name).suffix.lower()
            dest_dir = destinations.get(ext)
            if dest_dir is None:
                continue  # ignore other file types

            target = dest_dir / entry.name
            if target.exists():
                continue  # don't overwrite, makes the function safe to re-run

            transfer(entry.path, target)
            counts[ext] += 1

    return counts