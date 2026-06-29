from colorama import Fore, Back, Style, just_fix_windows_console

just_fix_windows_console()

def header(message: str):
    print(Style.BRIGHT + Fore.GREEN + message + Style.RESET_ALL)
    
def error(message: str):
    print(Style.BRIGHT + Fore.RED + message + Style.RESET_ALL)