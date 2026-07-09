import argparse
from pathlib import Path

import app as app
import jscb as jscb
import db as db

def main(args):
    if args.command == "db":
        print(db)
    if args.command == "set_api_key":
        data = {}
        with open(".env", 'r') as env:
            for line in env:
                data[line.split("=")[0].strip()] = line.split("=")[1].strip()
        
        data["NCBI_API_KEY"] = args.key
        
        with open(".env", 'w') as env:
            for key, value in data.items():
                env.write(f"{key}={value}\n")
                
    if args.command == "set_email":
        data = {}
        with open(".env", 'r') as env:
            for line in env:
                data[line.split("=")[0].strip()] = line.split("=")[1].strip()
        
        data["email"] = args.email
        
        with open(".env", 'w') as env:
            for key, value in data.items():
                env.write(f"{key}={value}\n")
                
    if args.command == "app":
        app.run(
            args.accession,
            args.db_path,
            args.dir,
            args.config,
            "diamond" if args.DIAMOND else "blastp",
            args.threads,
            args.sid,
            args.min_db,
            args.mge,
            args.fasta,
            args.add_new_to_db
        )
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers(dest="command", required=True)
    
    db_parser = subparser.add_parser("db", help="create/manage database")
    
    app_parser = subparser.add_parser("app", help="Run Alienness by Phyletic Pattern (APP) Algorithm")
    app_parser.add_argument("accession", help = "RefSeq/GenBank Accession Number")
    app_parser.add_argument("db_path", help = "Database Path", type=Path)
    app_parser.add_argument("-d", "--dir", help = "Working dir, default current", type=Path, default="/")
    app_parser.add_argument("-D", "--DIAMOND", help = "Use DIAMOND instead of BLAST", action="store_true")
    app_parser.add_argument("-t", "--threads", help = "Number of threads, default 1", type=int, default = 1)
    app_parser.add_argument('--sid', help = "Species ID override")
    app_parser.add_argument("-M", "--min_db", help = "Use smaller database, based off of original APP's database", action="store_true")
    app_parser.add_argument("-m", "--mge", help = "Do Marker Gene Enrichment", action="store_true")
    app_parser.add_argument("-f", "--fasta", help = "Query Fasta File", type=Path)
    app_parser.add_argument("-a", "--add_new_to_db", help = "Add any downloaded genomes to permanent database", action="store_true")  
    app_parser.add_argument("-c", "--config", help = "APP Config file. If not provided, will search cwd", type=Path, default="config.yaml")
    
    jscb_parser = subparser.add_parser("jscb", help="Run Jenson-Shannon Codon Bias (JSCB) Algorithm")
    
    set_api_key_parser = subparser.add_parser("set_api_key", help="Set NCBI API Key")
    set_api_key_parser.add_argument("key", help="NCBI API Key")
    
    set_email_parser = subparser.add_parser("set_email", help="Set email for NCBI reporting")
    set_email_parser.add_argument("email", help="Email")
    
    args = parser.parse_args()
    
    main(args)
