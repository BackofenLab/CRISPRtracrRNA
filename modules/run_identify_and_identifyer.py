import os
import subprocess
from os import listdir
from os.path import isfile, join
import shutil


def run_crispr_identify(folder_path, result_folder_path):
    """
    Runs the CRISPR identify on the given folder.
    """

    cmd = f"tools/CRISPRcasIdentifier/CRISPRcasIdentifier/crispr_identify.py --input_folder {folder_path} --result_folder {result_folder_path} --fast_run True"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()

def run_crispr_cas_identifier(folder_path, result_folder_path):
    """
    Runs the CRISPR cas identifyer on the given folder.
    """

    cmd = f"tools/CRISPRidentify/CRISPRidentify.py --input_folder {folder_path} --result_folder {result_folder_path} --fast_run True"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()


#          Components cas genes
######################################################
######################################################

def parse_csv_protein_file(file_name):
    dict_file_csv = {}
    with open(file_name) as f:
        lines = f.readlines()[1:]
    for line in lines:
        line_info = line.split(",")
        start = int(line_info[1])
        end = int(line_info[2])
        annotation = line_info[5].strip()
        key = (start, end)
        if "cas" in annotation:
            dict_file_csv[key] = annotation

    return dict_file_csv


def cas_identifier_result_folder_parser(folder_path):
    dict_cas_proteins = {}
    onlyfiles = [f for f in listdir(folder_path) if isfile(join(folder_path, f))]
    protein_files = [f for f in onlyfiles if "annotated_proteins" in f]
    for file_name in protein_files:
        full_path = join(folder_path, file_name)
        dict_file = parse_csv_protein_file(full_path)
        dict_cas_proteins = {**dict_cas_proteins, **dict_file}
    return dict_cas_proteins


def run_cas_identifier(file_name):
    try:
        cmd1 = "mkdir output_cas"
        cmd2 = "mkdir output_cas/cassette"
        process = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        process = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()
    except Exception:
        pass

    try:
        cmd = "mkdir output"
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()
    except Exception:
        pass

    command = f"python tools/CRISPRcasIdentifier/CRISPRcasIdentifier/CRISPRcasIdentifier.py -f {file_name} -ho output_cas/hmmsearch -st dna -co output_cas/cassette"
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    a, b = process.communicate()


def complete_info_with_cas_identifier(file_name):
    run_cas_identifier(file_name)
    dict_cas = cas_identifier_result_folder_parser("output_cas/cassette")
    try:
        shutil.rmtree("output_cas")
    except Exception:
        pass

    try:
        shutil.rmtree("output")
    except Exception:
        pass

    folder = "/".join(file_name.split("/")[:-1])
    file_base = file_name.split("/")[-1].split(".")[0]
    log_prodigal_file = file_base + "_prodigal.log"
    prodigal_prots = file_base + "_proteins.fa"
    full_path_log = join(folder, log_prodigal_file)
    full_path_prots = join(folder, prodigal_prots)

    try:
        os.remove(full_path_log)
    except Exception:
        pass

    try:
        os.remove(full_path_prots)
    except Exception:
        pass

    return dict_cas

