from colorama import Fore, Back, Style, just_fix_windows_console
from pathlib import Path
import yaml

just_fix_windows_console()

def header(message: str):
    print(Style.BRIGHT + Fore.GREEN + message + Style.RESET_ALL)
    
def error(message: str):
    print(Style.BRIGHT + Fore.RED + message + Style.RESET_ALL)

def update_yaml(path: Path, updates: dict) -> None:
    # Load existing contents (empty file returns None)
    data = {}
    if path.exists():
        with open(path) as f:
            data = yaml.safe_load(f) or {}

    # Modify
    data.update(updates)

    # Write back
    with open(path, "w") as f:
        yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)