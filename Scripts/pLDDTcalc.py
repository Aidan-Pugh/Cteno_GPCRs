### THIS SCRIPT ONLY WORKS FOR ALPHA_FOLD PROTEINS, NOT EXPERIMENTALLY PRODUCED PROTEINS

from collections import Counter
import os
import pandas as pd
import logging
import argparse
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

def configure_logging(unique_factor, header):
    """
    Function that sets up the log file so that everything the script does can be written down
    Requires the user to add a name for the log file
    
    """
    name_of_log_file = f"pLDDTcalc{header}{unique_factor}.log"

    logging.basicConfig(
        filename=name_of_log_file,
        level=logging.INFO,
        format ="%(message)s"
    )

    logging.info(f"{name_of_log_file}")
    logging.info("<>-------------<>-------------<>---------------<>--------------<>--------------<>------------<>-----------<>")
    print(f"LOGGING: {name_of_log_file} configured, information will be written here\n")
    return name_of_log_file

def get_plddt(pdb_file_name, directory, pre_cursor):
    """
    calculate the proportion of different plddt values
    """
    if not pdb_file_name or not pdb_file_name.endswith(".pdb"):
        # If the file does not exist, return an empty row (or a default row)
        return None
    
    logging.info(f"Collecting quality score for {pdb_file_name}")
    file_path = os.path.join(directory, pdb_file_name)
    file_data = []
    plddt_cats = {
        "very_low": lambda x: x < 50,
        "low": lambda x: 50 <= x < 70,
        "high": lambda x: 70 <= x < 90,
        "very_high": lambda x: x >= 90
    }
    order = ["very_low", "low", "high", "very_high"]

    # Collect the scores in one place
    all_plddt_scores = []
    with open(file_path, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                plddt_value = line[60:66].strip()
                all_plddt_scores.append(plddt_value)
    counts = Counter({key: sum(1 for s in all_plddt_scores if func(float(s))) for key, func in plddt_cats.items()})
    ordered_counts = {key:counts.get(key, 0) for key in order}
    logging.info(f"    {ordered_counts}")

    # Calculations: percentages and average plddt
    percentages = {key: (value / len(all_plddt_scores)) * 100 for key, value in ordered_counts.items()}
    logging.info(f"    Percentages: {percentages}")
    average_plddt = average_plddt = sum([float(s) for s in all_plddt_scores]) / len(all_plddt_scores)
    logging.info(f"    Average pLDDT score: {average_plddt}")

    row = {
        "entry": f"{pre_cursor}{pdb_file_name}",
        "file_pdb": pdb_file_name,
        "very_low": percentages.get("very_low", 0),
        "low": percentages.get("low", 0),
        "high": percentages.get("high", 0),
        "very_high": percentages.get("very_high", 0),
        "average_plddt": average_plddt
    }
    file_data.append(row)

    return file_data

def locate_directory(directory):
    """
    Checks the location of the directory and curates a list of files
    """
    all_files = []

    for f in os.listdir(directory):
        file_path = os.path.join(directory, f)
        if os.path.isfile(file_path):
            all_files.append(f)
        
    logging.info(f"Files found in {directory}: {len(all_files)}")
    print(f"{len(all_files)} files found, curating pLDDT scores")
    return all_files

def main_script(directory, header):
    """
    Compile the functions into one script
    """
    # Set up pLDDT.csv
    unique_factor = datetime.now().strftime("%d%b%Y_%H%M")
    columns = ["entry","file_pdb", "very_low", "low", "high", "very_high", "average_plddt"]
    data_frame_plddt = pd.DataFrame(columns=columns)
    log_file = configure_logging(unique_factor, header)
    count = 0

    files_to_check=locate_directory(directory)
    for f in files_to_check:
        row=get_plddt(f, directory, pre_cursor=header)
        if row != None:
            data_frame_row = pd.DataFrame(row)
            data_frame_plddt = pd.concat([data_frame_plddt, data_frame_row], ignore_index=True)
            count = count + 1
            print(f"{count}/{len(files_to_check)}: {f}")

        if row == None:
            print(f"     X {f} has no available structure")
            logging.info(f"    No structure available for {f}")
            count = count + 1

    print(f"\n\nTo see unfound structures, search 'XX:' in {log_file}")
    data_frame_plddt.to_csv(f"plddtcalc{header}{unique_factor}.csv")




def main():
    reset = "\033[0m"
    magenta = "\033[95m"
    red = "\033[91m"
    ascii_banner = f"""
{magenta}<>{red}-----------{magenta}<>{red}-----------{magenta}<>{red}-----------{magenta}<>{red}-----------{magenta}<>
{red}   {magenta}_______{red}_{magenta}_____{red}____  ____ _____{magenta}_{red}________{magenta}_{red}______      
{red}    ____ | |   |  _ \\|  _ \\_   _|{magenta}__ __ _| | ___ 
{red}   |  _ \\| |   | | | | | | || |{magenta}/ __/ _` | |/ __|
{red}   | |_) | |___| |_| | |_| || |{magenta} (_| (_| | | (__ 
{red}   | .__/|_____|____/|____/ |_|{magenta}\\___\\__,_|_|\\___|
{red}   |_|{magenta}_________________________{red}_________________

{magenta}<>{red}-----------{magenta}<>{red}-----------{magenta}<>{red}-----------{magenta}<>{red}-----------{magenta}<>{reset}                                                        
"""
    usage_format = "Usage: python pLDDTcalc.py -e entries.txt"

    parser = argparse.ArgumentParser(
        description=(f"{ascii_banner}\n\nRetrieve 3D protein structure files, prioritising files from PDBe over AlphaFold, from UniProt IDs\n{usage_format}"),
        formatter_class=argparse.RawTextHelpFormatter,
        usage = argparse.SUPPRESS,)
    
    # EXPECTED ARGUMENTS:
    parser.add_argument("-d", type=str, required=True, 
                        help="Directory containing 3D structures")
    parser.add_argument("-f", type=str, required=False, default = "", 
                        help = "If the file names need anything appended to the start (e.g 7tm2_) then write it here")
    args = parser.parse_args()
    # Use the arguments to state how they will look in the script
    source_directory = args.d
    header = args.f
    
    # Your code logic here, preferably in the form of a whole-script function

    print(ascii_banner)
    main_script(source_directory, header)
    

if __name__ == "__main__":
    main()
