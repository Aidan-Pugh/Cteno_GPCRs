## SCRIPT TO RETRIEVE THE PDB files from uniprot or PDBe

# If the protein structure exists on PDBe, prioritze retrieving it from there,
# as this is likely experimentally derived. If not, retrieve AlphaFold
import pandas as pd
import os
import requests
import json
import argparse
from datetime import datetime
import logging
from collections import Counter
import warnings
warnings.filterwarnings("ignore")

def configure_logging(unique_factor):
    """
    Function that sets up the log file so that everything the script does can be written down
    Requires the user to add a name for the log file
    
    """
    name_of_log_file = "Log_" + unique_factor + ".log"

    logging.basicConfig(
        filename=name_of_log_file,
        level=logging.INFO,
        format ="%(message)s"
    )

    logging.info(f"{name_of_log_file}")
    logging.info("<>-------------<>-------------<>---------------<>--------------<>--------------<>------------<>-----------<>")
    print(f"LOGGING: {name_of_log_file} configured, information will be written here")
    return name_of_log_file

def read_entries(entry_txt):
    """
    Read in entries from a text file
    """
    print(f"READING: Creating list of entries from file {entry_txt}")

    list_of_entries = []
    with open(entry_txt, "r") as file:
        for line in file:
            list_of_entries.append(line.strip())

    logging.info(f"Number of files to retrieve: {len(list_of_entries)}")
    return list_of_entries

def create_directory(unique_factor):
    """
    Creates a directory where all the retrieved files will be stored.
    Name is unique, by using the date and time of collection to prevent overwriting.
    """
    new_directory = "PDB_collection_" + unique_factor
    os.makedirs(new_directory, exist_ok=True)
    logging.info(f"Files will be stored in {new_directory}")
    print(f"STORAGE: Downloaded files will be stored in {new_directory}")

    return new_directory

def get_pdb(entry):
    """
    Get the pdb from UniProt, assuming there is a reviewed/experimentally derived structure
    
    """
    url = f"https://rest.uniprot.org/uniprotkb/{entry}.json"
    response = requests.get(url)
    if response.status_code == 200: #This is a success
        data=response.json()
        # rest_info=json.dumps(data, indent=2)
        pdb_entries = data.get("uniProtKBCrossReferences", [])
        for db in pdb_entries:
            if db.get("database") == "PDBsum":
                pdb_id = db.get("id")
                logging.info(f"PD: {entry} has a non-AlphaFold structure available at ID: {pdb_id}")

                return pdb_id
    return None

def check_uniprot_ID(entry):
    """
    Run this is the get_pdb function fails to return an ID
    """
    url = f"https://rest.uniprot.org/uniprotkb/{entry}.json"
    response = requests.get(url)
    if response.status_code == 200: #This is a success
        data=response.json()

        AF_entries = data.get("uniProtKBCrossReferences", [])
        for db in AF_entries:
            if db.get("database") == "AlphaFoldDB":
                AF_id = db.get("id")

                logging.info(f"AF: {entry} has an AlphaFold structure available at ID: {AF_id}")
                return AF_id
            
    return None

def collect_ID(entry):
    """
    Combines to functions and returns the ID
    """
    # ID from pdb 
    ID_collected = get_pdb(entry)
    if ID_collected != None:
        entry_id = [entry, ID_collected]
        return entry_id
    # AlphaFold ID
    if ID_collected == None:
        ID_collected = check_uniprot_ID(entry)
        if ID_collected != None:
            entry_id = [entry, ID_collected]
            return entry_id
    # No 3D structure found
    if ID_collected == None:
        logging.info(f"XX: {entry} has no ID")
        entry_id = [entry, None]
        return entry_id
        
def download(entry_id, storage_location):
    """
    Download the files to a specified location
    """
    if entry_id[1] == None:
        return None
    # Non-AF strutures
    if len(entry_id[1]) == 4:
        url = f"https://ebi.ac.uk/pdbe/entry-files/download/pdb{entry_id[1].lower()}.ent"
        source = "PDB"
    else:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{entry_id[1]}-F1-model_v4.pdb"
        source = "AF"

    response = requests.get(url)
    if response.status_code == 200:
        file_path = os.path.join(storage_location, f"{entry_id[0]}_{source}.pdb")
        with open(file_path, "wb") as file:
            file.write(response.content)

    logging.info(f"    | {entry_id[0]}/{entry_id[1]} collected from {source}")
    logging.info(f"    | {entry_id[0]}_{source}.pdb")
    # Return the name of the file
    return f"{entry_id[0]}_{source}.pdb"

def get_plddt(pdb_file_name, directory):
    if pdb_file_name == None:
        # If the file doesn"t exist, return an empty row (or a default row)
        return None
    

    logging.info(f"    Collecting quality score for {pdb_file_name}")
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
    with open(file_path, "r") as pdb_files:
        for line in pdb_files:
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

    row = {
        "entry": pdb_file_name.split("_")[0],
        "file_pdb": pdb_file_name,
        "very_low%": percentages.get("very_low", 0),
        "low%": percentages.get("low", 0),
        "high%": percentages.get("high", 0),
        "very_high%": percentages.get("very_high", 0),
        "average_plddt": average_plddt
    }

    file_data.append(row)


    return file_data

def main_script(entries_txt):
    """
    Compile all the other functions into this script:
    """
    print("\n<>-------------<>--------------<> PDB COLLECTOR <>--------------<>-------------<>")
    # Set up save location
    unique_factor=datetime.now().strftime("%d%b%Y_%H%M")
    log_name=configure_logging(unique_factor)
    storage_location=create_directory(unique_factor)
    entries = read_entries(entries_txt)
    logging.info("<>-------------<>-------------<>---------------<>--------------<>--------------<>------------<>-----------<>")
    logging.info("Collecting files for the following accessions:\n")
    print("\nCollecting files and calculating plddt scores...")
    columns = ["entry","file_pdb", "very_low%", "low%", "high%", "very_high%", "average_plddt"]
    data_frame_plddt = pd.DataFrame(columns=columns)

    count = 0
    num_entries=len(entries)
    for e in entries:
        count = count + 1
        reference_entry_id = collect_ID(e)
        new_file_name=download(reference_entry_id, storage_location)
        print(f"{count}/{num_entries}: {new_file_name}")

        file_data_retrieved=get_plddt(new_file_name, storage_location)
        if file_data_retrieved != None:
            data_frame_row = pd.DataFrame(file_data_retrieved)
            data_frame_plddt = pd.concat([data_frame_plddt, data_frame_row], ignore_index=True)
        if file_data_retrieved == None:
            print(f"     X {e} has no available structure")
            logging.info(f"    No structure available for {e}")
    
    print(f"\n\nTo see unfound structures, search 'XX:' in {log_name}")
    data_frame_plddt.to_csv(f"0_plddt_scores_{storage_location}.csv")

    logging.info("\n<>-------------<>--------------<>{ PDB files collected }<>-------------<>--------------<>")
    print("\nPLDDT: Percentages of each bin and average plddt score calculated")
    return print("\n<>-------------<>--------------<> EVERYTHING HAS BEEN COLLECTED <>--------------<>-------------<>")



def main():
    cyan = "\033[96m"
    blue = "\033[94m"
    reset = "\033[0m"

    ascii_banner = f"""
{cyan}<>{blue}------------{cyan}<>{blue}------------{cyan}<>{blue}------------{cyan}<>{blue}------------{cyan}<>
 {blue}|  _ \\|  _ \\| __ ){cyan}  ___ ___ | | | ___  ___| |_ ___  _ __  
 {blue}| |_) | | | |  _ \\{cyan} / __/ _ \\| | |/ _ \\/ __| __/ _ \\| '__| 
 {blue}|  __/| |_| | |_) {cyan}| (_| (_) | | |  __/ (__| || (_) | |    
 {blue}|_|   |____/|____/{cyan} \\___\\___/|_|_|\\___|\\___|\\__\\___/|_|    
{cyan}<>{blue}------------{cyan}<>{blue}------------{cyan}<>{blue}------------{cyan}<>{blue}------------{cyan}<>\n{reset}                                                          
"""
    usage_format = "Usage: python PDBcollector.py -e entries.txt"

    parser = argparse.ArgumentParser(
        description=(f"{ascii_banner}\n\nRetrieve 3D protein structure files, prioritising files from PDBe over AlphaFold, from UniProt IDs\n{usage_format}"),
        formatter_class=argparse.RawTextHelpFormatter,
        usage = argparse.SUPPRESS,)
    
    # EXPECTED ARGUMENTS:
    parser.add_argument("-e", type=str, required=True, 
                        help="Entries list, as a text file. Each new line contains a UniProt ID")
    args = parser.parse_args()
    # Use the arguments to state how they will look in the script
    entry_txt = args.e
    # Your code logic here, preferably in the form of a whole-script function

    print(ascii_banner)
    main_script(entry_txt)
    

if __name__ == "__main__":
    main()
