#!/usr/bin/python3
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

Contact: Ethan Coppenbarger - ecoppen@iastate.edu

Usage:
$ python3.6 ScanFold-Fold2.0.py fasta_filename [options]

"""

import argparse
import logging
import numpy as np
import os
import re
import shutil
import sys
import time

from ScanFoldFunctions import *


# import RNAstructure


def main(args):
    start_time = time.time()
    filename = args.filename
    filter_value = int(args.f)
    competition = int(args.c)
    name = args.id
    global_refold = args.global_refold
    folder_name = args.folder_name
    temperature = int(args.t)
    algo = "rnafold"
    type = "mono"
    cwd = os.getcwd()
    fasta_file_path = args.fasta_file_path
    es_path = args.es_path
    igv_path_in = args.igv_path
    inforna_path_in = args.inforna_path
    
    # parse string output of glob.glob into list
    if filename[0] == '[':
        filename = eval(filename)

    # This sets strand to forward (but can be set up as an argument if needed)
    strand = 1
    for tsv_file in filename:
        os.chdir(cwd)
        with open(tsv_file, 'r') as f:
            # Read all lines from ScanFold-Scan file (including header)
            try:
                filename_split = tsv_file.split(".tsv")
                outname = filename_split[0]
            except:
                logging.inro("Input name should have .tsv extension")
                outname = tsv_file

            # handle case of input file being an absolute path
            outname = os.path.basename(outname)
            curr_dir = cwd
            if len(filename) > 1:
                unique_folder_name = str(outname + "_fold")
                curr_dir = os.path.join(cwd, unique_folder_name)
                if not os.path.exists(curr_dir):
                    os.mkdir(curr_dir)
                os.chdir(curr_dir)
            # create a folder for extracted_structures:
            # es_path = "extracted_structures"
            if not os.path.isabs(es_path):
                extract_path = os.path.join(curr_dir, es_path)
                if not os.path.exists(extract_path):
                    os.mkdir(extract_path)

            # create a folder for igv files:
            # igv_path = "igv_files"
            if not os.path.isabs(igv_path_in):
                igv_path = os.path.join(curr_dir, igv_path_in)
                if not os.path.exists(igv_path):
                    os.mkdir(igv_path)

            # create a folder for inforna_structures:
            # inforna_directory = "inforna_structures"
            if not os.path.isabs(inforna_path_in):
                inforna_path = os.path.join(curr_dir, inforna_path_in)
                if not os.path.exists(inforna_path):
                    os.mkdir(inforna_path)
            if folder_name is not None:
                try:
                    if len(filename) > 1:
                        logging.info("Putting output in folder named:" + unique_folder_name)
                        name = unique_folder_name
                        full_output_path = os.path.join(curr_dir, unique_folder_name)
                    else:
                        logging.info("Putting output in folder named:" + folder_name)
                        name = folder_name
                        full_output_path = os.path.join(curr_dir, folder_name)
                    os.mkdir(full_output_path)
                    os.chdir(full_output_path)
                    folder_name = str(folder_name)
                except:
                    logging.info(f"WARNING: FOLDER NAME [{folder_name}] IS NOT UNIQUE")
                    raise FileExistsError(
                        "Folder name exists. Manually input unique folder name via --folder_name [FOLDER NAME]")
            else:
                try:
                    folder_name = outname + ".Fold_output"
                    if len(filename) > 1:
                        logging.info("Putting output in folder named:" + unique_folder_name)
                        name = unique_folder_name
                        full_output_path = os.path.join(curr_dir, unique_folder_name)
                    else:
                        logging.info("Putting output in folder named:" + folder_name)
                        full_output_path = os.path.join(curr_dir, folder_name)
                    os.mkdir(full_output_path)
                    os.chdir(full_output_path)
                    folder_name = str(outname + "_fold")
                except:
                    logging.info(f"WARNING: FOLDER NAME [{folder_name}] IS NOT UNIQUE")
                    raise FileExistsError(
                        "Folder name exists. Manually input unique folder name via --folder_name [FOLDER NAME]")
            lines = f.readlines()[1:]  # skip header
            steplines1 = int(lines[0].split("\t")[0])
            steplines2 = int(lines[1].split("\t")[0])
            step_size = steplines2 - steplines1
            # Initialize bp dictionary and z-score lists
            z_score_list = []
            bp_dict = {}

            # Reset file names #
            dbn_file_path = outname + ".AllDBN_refold.txt"
            dbn_file_path1 = outname + ".NoFilter"
            dbn_file_path2 = outname + ".Zavg_-1"
            dbn_file_path3 = outname + ".Zavg_-2"
            dbn_file_path4 = outname + ".AllDBN.txt"
            structure_extract_file = outname + ".structure_extracts.txt"


            # Generate nucleotide dictionary to assign each nucleotide in sequence a key
            nuc_dict = NucleotideDictionary(lines)
            full_fasta_sequence = nuc_dict_to_seq(nuc_dict)
            #print("original full_fasta_sequence: " + str(full_fasta_sequence))
            logging.info("Sequence length: " + str(len(nuc_dict)) + "nt")
            # if len(nuc_dict) > 1000:
            #     raise SystemExit('Input sequence is longer than 1000 nt (check: https://github.com/moss-lab/ScanFold)')

            # Determine start and end coordinate values
            start_coordinate = str(list(nuc_dict.keys())[0])
            # logging.info(start_coordinate)
            end_coordinate = str(list(nuc_dict.keys())[-1])
            # logging.info(end_coordinate)

            # initiate value lists
            mfe_list = []
            zscore_list = []
            pvalue_list = []
            ed_list = []
            # Iterate through input file, read each rows metrics, sequence, etc.
            logging.info("Reading sequence and structures...")
            for row in lines:

                # Ignore blank lines
                if not row.strip():
                    continue

                # Main loop to find all i-j pairs per i-nucleotide
                else:
                    data = row.split('\t')
                    icoordinate = data[0]
                    jcoordinate = data[1]
                    temp = data[2]
                    mfe = float(data[3])
                    zscore = float(data[4])
                    pvalue = data[5]
                    ed = float(data[6])
                    sequence_raw = simple_transcribe(str(data[7]))
                    structure_raw = str(data[8])
                    #print("Raw structure: " + str(structure_raw))

                    strand = 1

                    sequence = list(sequence_raw)
                    structure = list(structure_raw)

                    # Define window coordinates as string
                    # window = str(str(icoordinate)+"-"+str(jcoordinate))

                    # Determine length of window
                    length = len(sequence)

                    # Append window z-score to list (to calculate overall z-score)
                    z_score_list.append(zscore)

                    # Iterate through dot bracket structure to determine locations of bps
                    i = 0
                    while i < length:
                        # Unpaired nucleotide
                        if structure[i] == '.':
                            nucleotide = sequence[i]
                            coordinate = (i + int(icoordinate))
                            x = NucPair(nucleotide, coordinate, nucleotide, coordinate,
                                        zscore, mfe, ed)
                            try:
                                y = bp_dict[coordinate]
                                y.append(x)
                            except:
                                bp_dict[coordinate] = []
                                y = bp_dict[coordinate]
                                y.append(x)
                            i += 1
                        # Paired nucleotide
                        else:
                            i += 1

                    # Initiate base pair tabulation variables
                    bond_order = []
                    bond_count = 0
                    bond_order_bracket = []
                    bond_count_bracket = 0
                    bond_order_curly = []
                    bond_count_curly = 0
                    bond_order_carrot = []
                    bond_count_carrot = 0
                    length = len(structure)
                    # Iterate through sequence to assign nucleotides to structure type
                    m = 0
                    n = 0
                    curly = 0
                    carrot = 0
                    paren_types = ['.','(',')','{','}','[',']','<','>']
                    while m < length:
                        if structure[m] == '(':
                            bond_count += 1
                            bond_order.append(bond_count)
                            m += 1
                        elif structure[m] == ')':
                            bond_order.append(bond_count)
                            bond_count -= 1
                            m += 1
                        elif structure[m] in paren_types:
                            bond_order.append(0)
                            m += 1
                        else:
                            print(structure)

                            print("Error in asigning nucleodides () to structure type")
                            sys.exit("line 222")
                    while n < length:
                        if structure[n] == '[':
                            bond_count_bracket += 1
                            bond_order_bracket.append(bond_count_bracket)
                            n += 1
                        elif structure[n] == ']':
                            bond_order_bracket.append(bond_count_bracket)
                            bond_count_bracket -= 1
                            n += 1
                        elif structure[n] in paren_types:
                            bond_order_bracket.append(0)
                            n += 1
                        else:
                            print("Error in asigning nucleodides [] to structure type")
                            print(structure)
                            sys.exit()

                    while curly < length:
                        if structure[curly] == '{':
                            bond_count_curly += 1
                            bond_order_curly.append(bond_count_curly)
                            curly += 1
                        elif structure[curly] == '}':
                            bond_order_curly.append(bond_count_curly)
                            bond_count_curly -= 1
                            curly += 1
                        elif structure[curly] in paren_types:
                            bond_order_curly.append(0)
                            curly += 1
                        else:
                            print("Error in asigning nucleodides {} to structure type")
                            print(structure)
                            sys.exit()

                    while carrot < length:
                        if structure[carrot] == '<':
                            bond_count_carrot += 1
                            bond_order_carrot.append(bond_count_carrot)
                            carrot += 1
                        elif structure[carrot] == '>':
                            bond_order_carrot.append(bond_count_carrot)
                            bond_count_carrot -= 1
                            carrot += 1
                        elif structure[carrot] in paren_types:
                            bond_order_carrot.append(0)
                            carrot += 1
                        else:
                            print("Error in asigning nucleodides to structure type")
                            print(structure)
                            sys.exit()
                    # Initiate base_pair list
                    base_pairs = []

                    # Create empty variable named test
                    # test = ""

                    # Iterate through bond order
                    j = 0
                    while j < len(bond_order):
                        if bond_order[j] != 0:
                            test = bond_order[j]
                            base_pairs.append(j + 1)
                            bond_order[j] = 0
                            j += 1
                            k = 0
                            while k < len(bond_order):
                                if bond_order[k] == test:
                                    base_pairs.append(k + 1)
                                    bond_order[k] = 0
                                    k += 1
                                else:
                                    k += 1
                        elif bond_order_bracket[j] != 0:
                            test = bond_order_bracket[j]
                            base_pairs.append(j+1)
                            bond_order_bracket[j] = 0
                            j += 1
                            k = 0
                            while k < len(bond_order_bracket):
                                if bond_order_bracket[k] == test:
                                    base_pairs.append(k+1)
                                    bond_order_bracket[k] = 0
                                    k += 1
                                else:
                                    k += 1

                        elif bond_order_curly[j] != 0:
                            test = bond_order_curly[j]
                            base_pairs.append(j+1)
                            bond_order_curly[j] = 0
                            j += 1
                            k = 0
                            while k < len(bond_order_curly):
                                if bond_order_curly[k] == test:
                                    base_pairs.append(k+1)
                                    bond_order_curly[k] = 0
                                    k += 1
                                else:
                                    k += 1

                        elif bond_order_carrot[j] != 0:
                            test = bond_order_carrot[j]
                            base_pairs.append(j+1)
                            bond_order_carrot[j] = 0
                            j += 1
                            k = 0
                            while k < len(bond_order_carrot):
                                if bond_order_carrot[k] == test:
                                    base_pairs.append(k+1)
                                    bond_order_carrot[k] = 0
                                    k += 1
                                else:
                                    k += 1
                        else:
                            j += 1
                    # Iterate through "base_pairs" "to define bps
                    l = 0
                    while l < len(base_pairs):
                        lbp = base_pairs[l]
                        rbp = base_pairs[l + 1]

                        lb = str(sequence[int(lbp) - 1])
                        rb = str(sequence[int(rbp) - 1])

                        lbp_coord = int(int(lbp) + int(icoordinate) - 1)
                        rbp_coord = int(int(rbp) + int(icoordinate) - 1)
                        x = NucPair(lb, lbp_coord, rb, rbp_coord, zscore, mfe, ed)
                        z = NucPair(rb, rbp_coord, lb, lbp_coord, zscore, mfe, ed)

                        # Try to append i-j pair to i-nuc for left i-nuc
                        try:
                            y = bp_dict[lbp_coord]
                            y.append(x)
                        # If i-nuc not defined, define it
                        except:
                            bp_dict[lbp_coord] = []
                            y = bp_dict[lbp_coord]
                            y.append(x)

                        # Try to append i-j pair to i-nuc for right i-nuc
                        try:
                            w = bp_dict[rbp_coord]
                            w.append(z)
                        # If i-nuc not defined, define it
                        except:
                            bp_dict[rbp_coord] = []
                            w = bp_dict[rbp_coord]
                            w.append(z)
                        l += 2
                    mfe_list.append(mfe)
                    zscore_list.append(zscore)
                    pvalue_list.append(pvalue)
                    ed_list.append(ed)
                # Define OVERALL values of metrics
                meanz = float(np.mean(z_score_list))
                sdz = float(np.std(z_score_list))
                minz = min(z_score_list)
                stdz = float(np.std(z_score_list))

                one_sig_below = float(meanz - stdz)
                two_sig_below = float(meanz - (2 * stdz))

        # Initiate global dictionaries to store best base pairs
        best_bps = {}
        best_sum_bps = {}
        best_sum_bps_means = {}
        best_total_window_mean_bps = {}

        # Iterate through initial i-nuc dictionary to determine best base pairs (round 1)
        elapsed_time = round((time.time() - start_time), 2)
        logging.info("Elapsed time: " + str(elapsed_time) + "s")
        logging.info("Determining best base pairs...")
        for k, v in sorted(bp_dict.items()):
            # Initiate local dictionaries to store metrics per nucleotide
            zscore_dict = {}
            pair_dict = {}
            mfe_dict = {}
            ed_dict = {}

            # Iterate through all i-j pairs per i-nucleotide to store metrics for each
            for pair in v:
                # Create a key  which contains nucleotide and coordinate info
                partner_key = str(pair.jnucleotide) + "-" + str(pair.jcoordinate)
                # logging.info(partner_key)

                # Create a variable which contains all i-j pair info
                x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                            pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)

                # Try to append the value of each metric to metric lists per i-nuc
                try:
                    y = zscore_dict[partner_key]
                    y.append(pair.zscore)

                    m = mfe_dict[partner_key]
                    m.append(pair.mfe)

                    e = ed_dict[partner_key]
                    e.append(pair.ed)

                    z = pair_dict[partner_key]
                    z.append(x)

                # If pair not defined, define it
                except:
                    zscore_dict[partner_key] = []
                    y = zscore_dict[partner_key]
                    y.append(pair.zscore)
                    pair_dict[partner_key] = []

                    mfe_dict[partner_key] = []
                    m = mfe_dict[partner_key]
                    m.append(pair.mfe)

                    ed_dict[partner_key] = []
                    e = ed_dict[partner_key]
                    e.append(pair.ed)

                    z = pair_dict[partner_key]
                    z.append(x)

            # Calculate and store sum of z-score per i-j pair
            sum_z = {}
            sum_z_lengths = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = sum(v1)
                test = sum_z[k1] = sum(v1)
                sum_z_lengths[k1] = len(sum_z)

            # Calculate and store mean of z-score per i-j pair
            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = np.mean(v1)
                test = mean_z[k1] = np.mean(v1)

            # Calculate and store mean MFE per i-j pair
            mean_mfe = {}
            for k1, v1 in mfe_dict.items():
                mean_mfe[k1] = np.mean(v1)

            # Calculate and store mean ED per i-j pair
            mean_ed = {}
            for k1, v1 in ed_dict.items():
                mean_ed[k1] = np.mean(v1)

            # Calculate and store total window counts per i-j pair
            total_windows = 0
            num_bp = 0
            for k1, v1 in zscore_dict.items():
                total_windows = total_windows + len(v1)
                key_data = re.split("-", str(k1))
                key_i = str(key_data[1])
                if int(k) == int(key_i):
                    continue
                if int(k) != int(key_i):
                    num_bp += 1

            # Print first line of log file tables (first half of log file)
            k_nuc = str(nuc_dict[k].nucleotide)
            # logging.info(k_nuc)
            with open(outname + ".ScanFold.log", 'a+', newline='\n') as log_total:
                log_total.write("\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED\tSumZ\tSumZ/#TotalWindows\tBPs= " + str(num_bp) + "\n")
                log_total.write("nt-" + str(k) + "\t-\t" + str(k_nuc) + "\t" + str(total_windows) + "\t-\t-\t-\t-\t-" + "\n")
                # Print remainder of log file tables (first half of log file)
                total_window_mean_z = {}
                for k1, v1 in zscore_dict.items():
                    bp_window = str(len(v1))
                    key_data = re.split("-", str(k1))
                    key_nuc = str(key_data[0])
                    key_i = str(key_data[1])
                    total_window_mean_z[k1] = (sum(v1)) / total_windows
                    z_sum = str(round(sum(v1), 2))
                    z_avg = str(round(np.mean(v1), 2))
                    test = str(round(total_window_mean_z[k1], 2))
                    k1_mean_mfe = str(round(mean_mfe[k1], 2))
                    k1_mean_ed = str(round(mean_ed[k1], 2))
                    if int(k) == int(key_i):
                        # logging.info("iNuc is "+str(key_i))
                        log_total.write(str(k) + "\tNoBP\t" + key_nuc + "\t" + bp_window + "\t" + k1_mean_mfe + "\t" + z_avg + "\t" + k1_mean_ed + "\t" + z_sum + "\t" + str(test) + "\n")
                    else:
                        # logging.info("j is "+str(k))
                        log_total.write(str(k) + "\t" + key_i + "\t" + key_nuc + "\t" + bp_window + "\t" + k1_mean_mfe + "\t" + z_avg + "\t" + k1_mean_ed + "\t" + z_sum + "\t" + str(test) + "\n")
            # Define best_bp_key based on coverage-normalized z-score
            best_bp_key = min(total_window_mean_z, key=total_window_mean_z.get)

            # Access best i-j NucPairs for each metric using best_bp_key
            best_bp_mean_z = mean_z[best_bp_key]
            best_bp_sum_z = sum_z[best_bp_key]
            best_bp_mean_mfe = mean_mfe[best_bp_key]
            best_bp_mean_ed = mean_ed[best_bp_key]
            best_total_window_mean_z = total_window_mean_z[best_bp_key]

            #Access best i-j pair info from key name
            best_bp_data = re.split("-", best_bp_key)
            best_nucleotide = best_bp_data[0]
            best_coordinate = best_bp_data[1]

            #Fill dictionary with coverage normalized z-score
            #logging.info("Determining best base pair for nucleotide ", k)
            best_total_window_mean_bps[k] = (NucPair((nuc_dict[k]).nucleotide,
                                            nuc_dict[k].coordinate, best_nucleotide,
                                            best_coordinate, best_total_window_mean_z,
                                            best_bp_mean_mfe, best_bp_mean_ed))

            #Fill dictionary with coverage average z-score
            best_bps[k] = (NucPair((nuc_dict[k]).nucleotide, (nuc_dict[k]).coordinate,
                                    best_nucleotide, best_coordinate, best_bp_mean_z,
                                    best_bp_mean_mfe, best_bp_mean_ed))

        # Detect competing partners, and select final i-j pairs
        final_partners = {}
        elapsed_time = round((time.time() - start_time), 2)
        logging.info("Elapsed time: " + str(elapsed_time) + "s")

                # print header for final partner log file (log_win)
        with open(outname + ".ScanFold.FinalPartners.txt", 'a+', newline='\n') as log_win:
            log_win.write("\ni\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED\t*Indicates most favorable bp has competition; bp(j) has more favorable partner or is more likely to be unpaired\n")

        # Iterate through round 1 i-j pairs
        if competition == 1:
            # logging.info(start_coordinate, end_coordinate)
            logging.info("Detecting competing pairs...")
            j_coord_list = []
            # for k, v in sorted(best_bps.items()):
            #     logging.info(jcoordinate)
            #     j_coord_list.append(int(v.jcoordinate))

            for k, v in sorted(best_bps.items()):
                # logging.info(k, v.icoordinate, v.jcoordinate)
                test_k = int(k)
                # logging.info(sum(test_k == int(v.jcoordinate) for v in best_bps.values()))
                if sum(test_k == int(v.jcoordinate) for v in best_bps.values()) >= 0:
                    #print(start_coordinate, end_coordinate)
                    # Scan the entire dictionary:
                    # keys = range(int(start_coordinate), int(end_coordinate))

                    # Scan two window lengths flanking nucleotide:
                    #print('best bps: ' + str(len(best_bps)))
                    # print(length*4)
                    length = len(nuc_dict)
                    #print('length: ' + str(length))
                    if (len(best_bps) < length*4):
                        # print("Scanning full dictionary")
                        # Length of input less than length of flanks
                        # keys = range(int(start_coordinate), int(end_coordinate))
                        subdict = best_total_window_mean_bps

                    elif (
                        (v.icoordinate - length*(2)) >= int(start_coordinate) and
                        (v.icoordinate + (length*2)) <= int(end_coordinate)
                        ):
                        #print(str(v.icoordinate - length*(2)))
                        # print("MIDDLE")
                        keys = range(int(v.icoordinate-(length*2)), int(v.icoordinate+(length*2)))
                        subdict = {k: best_total_window_mean_bps[k] for k in keys}

                    elif (
                        int(v.icoordinate + (length*(2))) <= int(end_coordinate)and
                        (v.icoordinate + (length*2)) <= int(end_coordinate)
                    ):
                        #print("BEGINING"+str(v.icoordinate - (length*(2)))+" "+str(end_coordinate))
                        keys = range(int(start_coordinate), int(v.icoordinate+(length*2))+1)
                        subdict = {k: best_total_window_mean_bps[k] for k in keys}

                    elif (v.icoordinate + (length*2)) >= int(end_coordinate):
                        if v.icoordinate-(length*2) > 0:
                            #print("END"+str(v.icoordinate + (length*2)))
                            keys = range(int(v.icoordinate-(length*2)), int(end_coordinate))
                        else:
                            keys =range(int(v.icoordinate-(length*2)), int(end_coordinate))
                            subdict = {k: best_total_window_mean_bps[k] for k in keys}

                    elif len(best_bps) < length:
                            subdict = best_total_window_mean_bps

                    else:
                        print("Sub-dictionary error")
                        raise ValueError("Sub-dictionary error")

                    if len(subdict) >= 0:

                        # logging.info("Found competing pair for "+str(k))
                        # elapsed_time = round((time.time() - start_time), 2)
                        # logging.info(elapsed_time)
                        # logging.info("Detecting competing pairs for nuc ", k)
                        # For each i and j in i-j pair, detect competing pairs and append to dict
                        comp_pairs_i = competing_pairs(subdict, v.icoordinate)
                        # logging.info("i-pairs="+str(len(comp_pairs_i)))
                        comp_pairs_j = competing_pairs(subdict, v.jcoordinate)
                        # logging.info("j-pairs="+str(len(comp_pairs_j)))
                        total_pairs = []
                        # Put pairs competing with i from i-j pair into total pair dict for i-nuc
                        for key, pair in comp_pairs_i.items():
                            # logging.info("checking competing pairs for i")
                            # if k == 216:
                            #     logging.info(k, pair.icoordinate, pair.jcoordinate, pair.zscore)
                            total_pairs.append(competing_pairs(subdict,
                                                            pair.jcoordinate))
                            total_pairs.append(competing_pairs(subdict,
                                                            pair.icoordinate))
                        #Put pairs competing with j from i-j pair into total pair dict for i-nuc
                        for key, pair in comp_pairs_j.items():
                            # logging.info("checking competing pairs for j")
                            # if k == 216:
                            #     logging.info(k, pair.icoordinate, pair.jcoordinate, pair.zscore)
                            total_pairs.append(competing_pairs(subdict,
                                                            pair.jcoordinate))
                            total_pairs.append(competing_pairs(subdict,
                                                            pair.icoordinate))
                        # logging.info(str(k)+"nt Total comp pairs="+str(len(total_pairs)))

                        # Merge all dictionaries
                        merged_dict = {}
                        i = 0
                        for d in total_pairs:
                            # logging.info("merging competing dictionaries "+str(i))
                            for k1, v1 in d.items():
                                # if k == 216:
                                #     logging.info(k, k1, v1.icoordinate, v1.jcoordinate, v1.zscore)
                                merged_dict[i] = v1
                                i += 1

                        # #logging.info("MergedDict length for "+str(k)+"="+str(len(merged_dict)))
                        # #initiate best_basepair function, return best_bp based on sum
                        # if len(merged_dict) > 2:
                        #     bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                        #     #logging.info(str(len(merged_dict))+"__11111")
                        # else:
                        #     #logging.info("Nucleotide "+str(k))
                        #     bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                        #     #logging.info(str(len(merged_dict))+"____222222")
                        #     #bp = best_total_window_mean_bps[k]
                        # #Check if best basepair was connected to i-nucleotide (i.e., "k")
                        # logging.info(len(merged_dict))
                        if len(merged_dict) > 0:
                            bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                        else:
                            bp = NucPair(v.inucleotide,
                                         v.icoordinate,
                                         v.inucleotide,
                                         v.icoordinate,
                                         v.zscore,
                                         v.mfe,
                                         v.ed)

                        if (int(k) != bp.icoordinate) and (int(k) != int(bp.jcoordinate)):
                            # logging.info("1 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                            # if there was a competing i-j pair print it to log file instead:
                            with open(outname + ".ScanFold.FinalPartners.txt", 'a+', newline='\n') as log_win:
                                log_win.write("nt-" + str(k) + "*:\t" + str(v.icoordinate) + "\t" + v.jcoordinate + "\t" + str(round(v.mfe, 2)) + "\t" + str(round(v.zscore, 2)) + "\t" + str(round(v.ed, 2)) + "\n")
                            final_partners[k] = NucPair(v.inucleotide,
                                                        v.icoordinate,
                                                        v.inucleotide,
                                                        v.icoordinate,
                                                        best_bps[bp.icoordinate].zscore,
                                                        bp.mfe,
                                                        bp.ed)
                            #print(str(final_partners[k].inucleotide))
                            #print(str(final_partners[k].icoordinate))
                            #print(str(final_partners[k].jnucleotide))
                            #print(str(final_partners[k].jcoordinate))
                            #print(str(final_partners[k].zscore))
                            #print(str(final_partners[k].mfe))
                            #print(str(final_partners[k].ed))
                        #
                        # elif (int(v.icoordinate) == int(v.jcoordinate)) and (int(bp.icoordinate) != int(bp.jcoordinate)):
                        #     #Check for instance where competing base pair
                        #     logging.info(
                        #     "2 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+
                        #     " AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate)
                        #     )
                        #     log_win.write("nt-"+str(k)+"*:\t"+str(bp.icoordinate)+"\t"+bp.jcoordinate+"\t" +str(round(bp.mfe, 2)) +"\t"+str(round(bp.zscore, 2)) +"\t"+str(round(bp.ed, 2))+"\n")
                        #     final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                        #                                 bp.jnucleotide, bp.jcoordinate,
                        #                                 best_bps[bp.icoordinate].zscore,
                        #                                 best_bps[bp.icoordinate].mfe,
                        #                                 best_bps[bp.icoordinate].ed)
                        #
                        #
                        else:
                            # logging.info(
                            # "3 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+
                            # " AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate)
                            # )
                            with open(outname + ".ScanFold.FinalPartners.txt", 'a+', newline='\n') as log_win:
                                log_win.write("nt-" + str(k) + ":\t" + str(bp.icoordinate) + "\t" + str(bp.jcoordinate) + "\t" + str(round(best_bps[k].mfe, 2)) + "\t" + str(round(best_bps[k].zscore, 2)) + "\t" + str(round(best_bps[k].ed, 2)) + "\n")
                            final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                                                        bp.jnucleotide, bp.jcoordinate,
                                                        best_bps[bp.icoordinate].zscore,
                                                        best_bps[bp.icoordinate].mfe,
                                                        best_bps[bp.icoordinate].ed)
                            #print(str(final_partners[k].inucleotide))
                            #print(str(final_partners[k].icoordinate))
                            #print(str(final_partners[k].jnucleotide))
                            #print(str(final_partners[k].jcoordinate))
                            #print(str(final_partners[k].zscore))
                            #print(str(final_partners[k].mfe))
                            #print(str(final_partners[k].ed))
                    else:
                        continue
                else:
                    final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                                v.jnucleotide, v.jcoordinate,
                                                best_bps[k].zscore,
                                                best_bps[k].mfe,
                                                best_bps[k].ed)
                    #print(str(final_partners[k].inucleotide))
                    #print(str(final_partners[k].icoordinate))
                    #print(str(final_partners[k].jnucleotide))
                    #print(str(final_partners[k].jcoordinate))
                    #print(str(final_partners[k].zscore))
                    #print(str(final_partners[k].mfe))
                    #print(str(final_partners[k].ed))
                    # logging.info("No competing pair found for ", k)
                    continue

        # Write CT files
        if competition == 1:
            logging.info("Trying to write CT files with -c option")
            elapsed_time = str(round((time.time() - start_time), 2)) + "s"
            logging.info("Elapsed time: " + elapsed_time)
            logging.info("Writing CT files")

            write_ct(final_partners, str(dbn_file_path1) + ".ct", float(10), strand, name, start_coordinate)
            write_ct(final_partners, str(dbn_file_path2) + ".ct", float(-1), strand, name, start_coordinate)
            write_ct(final_partners, str(dbn_file_path3) + ".ct", float(-2), strand, name, start_coordinate)
            if filter_value != -2:
                write_ct(final_partners, str(user_filter_dbn) + ".ct", filter_value, strand, name, start_coordinate)
                makedbn(user_filter_dbn, "Zavg" + str(filter_value))
            makedbn(dbn_file_path1, "NoFilter")
            makedbn(dbn_file_path2, "Zavg_-1")
            makedbn(dbn_file_path3, "Zavg_-2")
            write_bp(final_partners, igv_path + "/" + outname + ".bp", start_coordinate, name, minz)
            write_wig_dict(final_partners, igv_path + "/" + outname + ".zavgs.wig", name, step_size, str("zscore"))
            write_wig_dict(final_partners, igv_path + "/" + outname + ".mfe_avgs.wig", name, step_size, str("mfe"))
            write_wig_dict(final_partners, igv_path + "/" + outname + ".ed_avgs.wig", name, step_size, str("ed"))
            write_bp(best_bps, igv_path + "/" + outname + ".ALL.bp", start_coordinate, name, minz)
            """
            write_bp(final_partners, outname + ".bp", start_coordinate, name, minz)
            if args.bp_track is not None:
                write_bp(final_partners, args.bp_track, start_coordinate, name, minz)

            write_wig_dict(final_partners, outname + ".zavgs.wig", name, step_size, str("zscore"))
            if args.final_partners_wig is not None:
                write_wig_dict(final_partners, args.final_partners_wig, name, step_size, str("zscore"))

            write_wig_dict(final_partners, outname + ".mfe_avgs.wig", name, step_size, str("mfe"))
            if args.mfe_wig_file_path is not None:
                write_wig_dict(final_partners, args.mfe_wig_file_path, name, step_size, str("mfe"))

            write_wig_dict(final_partners, outname + ".ed_avgs.wig", name, step_size, str("ed"))
            if args.ed_wig_file_path is not None:
                write_wig_dict(final_partners, args.ed_wig_file_path, name, step_size, str("ed"))

            write_bp(best_bps, outname + ".ALL.bp", start_coordinate, name, minz)
            """
        elif competition == 0:
            if args.bp_track is not None:
                write_bp(best_bps, args.bp_track, start_coordinate, name, minz)
        else:
            raise ValueError("Competition value not properly set")

        write_fasta(nuc_dict, outname + ".fa", name)
        if fasta_file_path is not None:
            write_fasta(nuc_dict, fasta_file_path, name)

        write_fai(nuc_dict, outname + ".fai", name)
        if args.fasta_index is not None:
            write_fai(nuc_dict, args.fasta_index, name)

        logging.info("ScanFold-Fold analysis complete! Refresh page to ensure proper loading of IGV")
        merge_files(str(dbn_file_path4), str(dbn_file_path1 + ".dbn"), str(dbn_file_path2 + ".dbn"),
                    str(dbn_file_path3 + ".dbn"))

        """ Begin the structure extract process """
        if global_refold:
            #logging.info("Attempting global refold")
            # create a separate DBN file
            with open(str(dbn_file_path), "w+", newline='\n') as dbn_log_file:
                try:
                    # fold the full fasta input as a fold compound (full_fc) using model params (md)
                    logging.info("Refolding full sequence using ScanFold results as constraints...")
                    elapsed_time = round((time.time() - start_time), 2)
                    logging.info("Elapsed time: " + str(elapsed_time) + "s")
                    md = RNA.md()
                    md.temperature = int(temperature)

                    # refold from -1 constraints
                    fc = RNA.fold_compound(str(full_fasta_sequence), md)
                    with open(str(dbn_file_path2 + ".dbn"), "r") as dbn_file_filter1:
                        lines = dbn_file_filter1.readlines()
                        filter1constraints = str(lines[2].rstrip())
                        #print("Filter 1 constraints: " + str(filter1constraints))
                    fc.hc_add_from_db(filter1constraints)
                    (refolded_filter1_structure, refolded_filter1_MFE) = fc.mfe()

                    # refold from -2 constraints
                    fc = RNA.fold_compound(str(full_fasta_sequence), md)
                    with open(str(dbn_file_path3 + ".dbn"), "r") as dbn_file_filter2:
                        lines = dbn_file_filter2.readlines()
                        filter2constraints = str(lines[2].rstrip())
                        #print("Filter 2 constraints: " + str(filter2constraints))
                    fc.hc_add_from_db(filter2constraints)
                    (refolded_filter2_structure, refolded_filter2_MFE) = fc.mfe()

                    # refold the global structure
                    full_fc = RNA.fold_compound(str(full_fasta_sequence), md)
                    (full_structure, full_MFE) = full_fc.mfe()

                    # write refolded structures to file
                    dbn_log_file.write(">" + str(name) + "\tGlobal Full MFE=" + str(full_MFE) + "\n" + str(full_fasta_sequence) + "\n" + str(full_structure) + "\n")
                    dbn_log_file.write(">" + str(name) + "\tRefolded with -1 constraints MFE=" + str(refolded_filter1_MFE) + "\n" + str(full_fasta_sequence) + "\n" + str(refolded_filter1_structure) + "\n")
                    dbn_log_file.write(">" + str(name) + "\tRefolded with -2 constraints MFE=" + str(refolded_filter2_MFE) + "\n" + str(full_fasta_sequence) + "\n" + str(refolded_filter2_structure) + "\n")
                except:
                    logging.info("Automatic refold for " + cur_record.name + " failed. Run manually")

                # assign filter structure for motif extraction
                if args.extract == 1:
                    full_filter_structure = filter1constraints
                    #print(full_filter_structure)
                elif args.extract == 2:
                    full_filter_structure = filter2constraints
                    #print(full_filter_structure)
                else:
                    raise ValueError("Constraint value error")

        if not global_refold:
            # skip refold, assign filter structure for motif extraction
            if args.extract == 1:
                with open(dbn_file_path2 + ".dbn", "r") as dbn_file_filter1:
                    lines = dbn_file_filter1.readlines()
                    full_filter_structure = str(lines[2].rstrip())
                    #print(full_filter_structure)
            elif args.extract == 2:
                with open(dbn_file_path3 + ".dbn", "r") as dbn_file_filter2:
                    lines = dbn_file_filter2.readlines()
                    full_filter_structure = str(lines[2].rstrip())
                    #print(full_filter_structure)
            else:
                raise ValueError("Constraint value error")

        """ Find nested "(..)" pairs """
        length = len(full_fasta_sequence)
        length_st = len(full_filter_structure)
        #print(str(length))
        #print(full_fasta_sequence)
        #print(str(length_st))
        #print(full_filter_structure)
        if length != length_st:
            raise ValueError("Length of sequence and structure do not match")
        bond_count = 0
        bond_order = []
        nuc_dict_refold = {}
        for n, i in enumerate(full_filter_structure):
            if i == '(':
                bond_count += 1
                bond_order.append(bond_count)
                nuc_dict_refold[n] = NucStructure(bond_count,
                                                  (n + 1),
                                                  full_fasta_sequence[n],
                                                  full_filter_structure[n])
            elif i == ')':
                bond_order.append(bond_count)
                bond_count -= 1
                nuc_dict_refold[n] = NucStructure(bond_count,
                                                  (n + 1),
                                                  full_fasta_sequence[n],
                                                  full_filter_structure[n])
            elif i == '.' or '<' or '>' or '{' or '}' or '[' or ']':
                bond_order.append(0)
                nuc_dict_refold[n] = NucStructure(bond_count,
                                                  (n + 1),
                                                  full_fasta_sequence[n],
                                                  full_filter_structure[n])
            else:
                raise ValueError("Unrecognized structure in constraint file")
            # print(str(nuc_dict_refold[n].bond_order))
            # print(str(nuc_dict_refold[n].nucleotide))
            # print(str(nuc_dict_refold[n].structure))
            # print(str(nuc_dict_refold[n].coordinate))

        """ Repeat the process looking for non-nested "<..>" pairs """
        bond_count_pk = 0
        bond_order_pk = []
        nuc_dict_pk = {}
        # Iterate through sequence to assign nucleotides to structure type
        for n, i in enumerate(full_filter_structure):
            if i == '<':
                bond_count_pk += 1
                bond_order_pk.append(bond_count)
                nuc_dict_pk[n] = NucStructure(bond_count_pk,
                                              (n + 1),
                                              full_fasta_sequence[n],
                                              full_filter_structure[n])
            elif i == '>':
                bond_order_pk.append(bond_count)
                bond_count_pk -= 1
                nuc_dict_pk[n] = NucStructure(bond_count_pk,
                                              (n + 1),
                                              full_fasta_sequence[n],
                                              full_filter_structure[n])
            elif i == '.' or '(' or ')' or '{' or '}' or '[' or ']':
                bond_order_pk.append(0)
                nuc_dict_pk[n] = NucStructure(bond_count_pk,
                                              (n + 1),
                                              full_fasta_sequence[n],
                                              full_filter_structure[n])
            else:
                raise ValueError("Unrecognized structure in constraint file")
            # print(str(nuc_dict_pk[n].bond_order))
            # print(str(nuc_dict_pk[n].nucleotide))
            # print(str(nuc_dict_pk[n].structure))
            # print(str(nuc_dict_pk[n].coordinate))

        """ Repeat the process looking for non-nested "[..]" pairs """
        bond_count_pk = 0
        bond_order_pk = []
        nuc_dict_pk = {}
        # Iterate through sequence to assign nucleotides to structure type
        for n, i in enumerate(full_filter_structure):
            if i == '[':
                bond_count_pk += 1
                bond_order_pk.append(bond_count)
                nuc_dict_pk[n] = NucStructure(bond_count_pk,
                                              (n + 1),
                                              full_fasta_sequence[n],
                                              full_filter_structure[n])
            elif i == ']':
                bond_order_pk.append(bond_count)
                bond_count_pk -= 1
                nuc_dict_pk[n] = NucStructure(bond_count_pk,
                                              (n + 1),
                                              full_fasta_sequence[n],
                                              full_filter_structure[n])
            elif i == '.' or '(' or ')' or '{' or '}' or '<' or ']':
                bond_order_pk.append(0)
                nuc_dict_pk[n] = NucStructure(bond_count_pk,
                                              (n + 1),
                                              full_fasta_sequence[n],
                                              full_filter_structure[n])
            else:
                raise ValueError("Unrecognized structure in constraint file")
        """ Identify structural motifs """
        structure_count = 0
        structure_start = []
        structure_end = []
        for j, i in enumerate(full_filter_structure):
            if (nuc_dict_refold[j].bond_order == 1) and (nuc_dict_refold[j].structure == '('):
                structure_count += 1
                # print(nuc_dict_refold[j].coordinate)
                # print(nuc_dict_refold[j].nucleotide)
                # print(nuc_dict_refold[j].structure)
                # print(structure_count)
                structure_start.append(NucStructure(structure_count,
                                                    nuc_dict_refold[j].coordinate,
                                                    nuc_dict_refold[j].nucleotide,
                                                    nuc_dict_refold[j].structure))
            elif (nuc_dict_refold[j].bond_order == 0) and (nuc_dict_refold[j].structure == ')'):
                # print(nuc_dict_refold[j].coordinate)
                # print(nuc_dict_refold[j].nucleotide)
                # print(nuc_dict_refold[j].structure)
                # print(structure_count)
                structure_end.append(NucStructure(structure_count,
                                                  nuc_dict_refold[j].coordinate,
                                                  nuc_dict_refold[j].nucleotide,
                                                  nuc_dict_refold[j].structure))
            else:
                continue

        """ Repeat for non-nested """
        structure_count_pk = 0
        structure_start_pk = []
        structure_end_pk = []
        for j, i in enumerate(full_filter_structure):
            if (nuc_dict_pk[j].bond_order == 1) and (nuc_dict_pk[j].structure == '<'):
                structure_count_pk += 1
                # print(nuc_dict_pk[j].coordinate)
                # print(nuc_dict_pk[j].nucleotide)
                # print(nuc_dict_pk[j].structure)
                # print(structure_count_pk)
                structure_start_pk.append(NucStructure(structure_count_pk,
                                                       nuc_dict_pk[j].coordinate,
                                                       nuc_dict_pk[j].nucleotide,
                                                       nuc_dict_pk[j].structure))
            elif (nuc_dict_pk[j].bond_order == 0) and (nuc_dict_pk[j].structure == '>'):
                # print(nuc_dict_pk[j].coordinate)
                # print(nuc_dict_pk[j].nucleotide)
                # print(nuc_dict_pk[j].structure)
                # print(structure_count_pk)
                structure_end_pk.append(NucStructure(structure_count_pk,
                                                     nuc_dict_pk[j].coordinate,
                                                     nuc_dict_pk[j].nucleotide,
                                                     nuc_dict_pk[j].structure))
            else:
                continue

        """Add extracted structure(s) to list"""
        logging.info("Extracting structures")
        motif_count = 0
        offset = 0  # set value to add unpaired nucleotides to 5' and 3' ends
        extracted_structure_list = []
        for motif_index, i in enumerate(structure_start):
            motif_count += 1
            motif_start_coord = int(structure_start[motif_index].coordinate - offset)
            if motif_start_coord < 0: motif_start_coord = 0
            motif_end_coord = int(structure_end[motif_index].coordinate + offset)
            if motif_end_coord > int(length): motif_end_coord = int(length)
            motif_sequence = ""
            motif_structure = ""
            for nt_index, value in nuc_dict_refold.items():
                if motif_start_coord <= nt_index + 1 <= motif_end_coord:
                    motif_sequence += str(value.nucleotide)
                    motif_structure += str(value.structure)
            extracted_structure_list.append(ExtractedStructure(motif_count,
                                                               motif_sequence,
                                                               motif_structure,
                                                               motif_start_coord,
                                                               motif_end_coord))
            #print(extracted_structure_list[motif_count - 1].structure_count)
            #print(extracted_structure_list[motif_count - 1].sequence)
            #print(extracted_structure_list[motif_count - 1].structure)
            #print(extracted_structure_list[motif_count - 1].i)
            #print(extracted_structure_list[motif_count - 1].j)
            #print()

        """ Repeat for non-nested """
        for motif_index_pk, i in enumerate(structure_start_pk):
            motif_count += 1
            motif_start_coord = int(structure_start_pk[motif_index_pk].coordinate - offset)
            if motif_start_coord < 0: motif_start_coord = 0
            motif_end_coord = int(structure_end_pk[motif_index_pk].coordinate + offset)
            if motif_end_coord > int(length): motif_end_coord = int(length)
            motif_sequence = ""
            motif_structure = ""
            for nt_index, value in nuc_dict_pk.items():
                if motif_start_coord <= nt_index + 1 <= motif_end_coord:
                    motif_sequence += str(value.nucleotide)
                    motif_structure += str(value.structure)
            extracted_structure_list.append(ExtractedStructure(motif_count,
                                                               motif_sequence,
                                                               motif_structure,
                                                               motif_start_coord,
                                                               motif_end_coord))
            #print(extracted_structure_list[motif_count - 1].structure_count)
            #print(extracted_structure_list[motif_count - 1].sequence)
            #print(extracted_structure_list[motif_count - 1].structure)
            #print(extracted_structure_list[motif_count - 1].i)
            #print(extracted_structure_list[motif_count - 1].j)
            #print()

        """ Obtain statistics for motifs """
        zscore_total = []
        pvalue_total = []
        mfe_total = []
        ed_total = []
        with open(structure_extract_file, "w", newline='\n') as extract_file:
            #print("extracted_structure_list length: " + str(len(extracted_structure_list)))
            #extracted_structure_list[0].describeMe()
            #for i in extracted_structure_list:
            #    i.describeMe()
            #for i in extracted_structure_list[:]:
            #    i.describeMe()
            for motif in extracted_structure_list:
                # print(motif)
                # print(motif.structure_count)
                # print(motif.sequence)
                # print(motif.structure)
                # print(motif.i)
                # print(motif.j)
                frag = motif.sequence
                fc = RNA.fold_compound(str(frag))  # creates "Fold Compound" object
                fc.hc_add_from_db(str(motif.structure))
                fc.pf()  # performs partition function calculations
                frag_q = (RNA.pf_fold(str(frag)))  # calculate partition function "fold" of fragment
                (MFE_structure, mfe_calc) = fc.mfe()  # calculate and define variables for mfe and structure
                mfe_calc = round(mfe_calc, 2)
                mfe_total.append(mfe_calc)
                (centroid, distance) = fc.centroid()  # calculate and define variables for centroid
                ed_calc = round(fc.mean_bp_distance(), 2)  # this calculates ED based on last calculated partition function
                ed_total.append(ed_calc)
                seqlist = []  # creates the list we will be filling with sequence fragments
                seqlist.append(frag)  # adds the native fragment to list
                scrambled_sequences = scramble(frag, 100, type)
                seqlist.extend(scrambled_sequences)
                energy_list = energies(seqlist, temperature, algo)
                try:
                    zscore = round(zscore_function(energy_list, 100), 2)
                except:
                    zscore = zscore_function(energy_list, 100)
                zscore_total.append(zscore)
                pvalue = round(pvalue_function(energy_list, 100), 2)
                pvalue_total.append(pvalue)
                accession = str(name)
                # ps_title = f"motif_{motif_num} coordinates {es.i} - {es.j}"
                es_dbn_path = f"{extract_path}/{name}_motif_{motif.structure_count}.dbn"
                with open(es_dbn_path, 'w', newline='\n') as es_dbn:
                    es_dbn.write(f">{name}_motif_{motif.structure_count}_coordinates:{motif.i}-{motif.j}_zscore={zscore}\n{frag}\n{MFE_structure}")
                dbn2ct(es_dbn_path)
                os.rename(f"{extract_path}/{name}_motif_{motif.structure_count}.ct",
                          f"{inforna_path}/{name}_motif_{motif.structure_count}.ct")

                # Create postscript files
                RNA.PS_rna_plot_a(frag, MFE_structure, extract_path + "/motif_" + str(motif.structure_count) + ".ps", '',
                                  '')

                # Set extracted structures up as GFF format
                gff_attributes = f'motif_{motif.structure_count};sequence={motif.sequence};structure={str(motif.structure)};refoldedMFE={str(MFE_structure)};MFE(kcal/mol)={str(mfe_calc)};z-score={str(zscore)};ED={str(ed_calc)}'
                extract_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str("."), str("RNA_sequence_secondary_structure"), str(int(motif.i + 1)), str(int(motif.j + 1)), str("."), str("."), str("."), gff_attributes))

        # except:
        #     logging.info("Structure Extract failed for "+folder_name+", must extract manually.")
        #     continue
        # logging.info("TEST")
        shutil.make_archive(inforna_path, 'zip', inforna_path)
        elapsed_time = round((time.time() - start_time), 2)
        logging.info("Total runtime: " + str(elapsed_time) + "s")
        logging.info("ScanFold-Fold analysis complete! Output found in folder named: " + full_output_path)

        if args.webserver:
            make_tar(args.webserver, full_output_path)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('filename', nargs="+",
                        help='input filename')
    parser.add_argument('-f', type=int, default=-2,
                        help='filter value')
    parser.add_argument('-c', type=int, default=1,
                        help='Competition (1 for disallow competition, 0 for allow; 1 by default)')
    parser.add_argument('--id', type=str, default="UserInput",
                        help='Name or ID of sequence being analyzed. Default "UserInput"')
    parser.add_argument('--global_refold', action='store_true', default=False,
                        help='Global refold option. Refold full sequence using Zavg <-1 and <-2 base pairs')
    parser.add_argument('--folder_name', type=str,
                        help='Name of output folder (defaults to header name or date/time)', default=None)
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')
    parser.add_argument('--extract', type=int, default='2',
                        help='Extract structures from minus 1 or minus 2 dbn file (2 or 1); Default = 2')
    parser.add_argument('--es_path', type=str, default="extracted_structures",
                        help='name of extracted structures file')
    parser.add_argument('--igv_path', type=str, default="igv_files",
                        help='name of IGV file')
    parser.add_argument('--inforna_path', type=str, default="inforna_structures",
                        help='name of inforna file')

    # webserver stuff
    parser.add_argument('--logfile', default=sys.stdout, type=argparse.FileType('w', encoding='UTF-8'),
                        help='Path to write log file to.')
    parser.add_argument('--loglevel', default="INFO", type=str,
                        help='Log level.')
    parser.add_argument('--webserver', type=str,
                        help='If provided, the output folder is compressed into a tar.gz file and written to the path specified by this parameter')
    parser.add_argument('--fasta_file_path', type=str,
                        help='fasta_file path')
    parser.add_argument('--fasta_index', type=str,
                        help='fasta index file path')
    parser.add_argument('--bp_track', type=str,
                        help='bp_track_file path')
    parser.add_argument('--ed_wig_file_path', type=str,
                        help='ed_wig_file_path')
    parser.add_argument('--mfe_wig_file_path', type=str,
                        help='mfe_wig_file_path')
    parser.add_argument('--pvalue_wig_file_path', type=str,
                        help='pvalue_wig_file_path')
    parser.add_argument('--zscore_wig_file_path', type=str,
                        help='zscore_wig_file_path')
    parser.add_argument('--final_partners_wig', type=str,
                        help='final partners wig file path')

    args = parser.parse_args()

    loglevel = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(loglevel, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logging.basicConfig(stream=args.logfile, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=loglevel)

    try:
        main(args)
    except Exception as e:

        if args.webserver:
            # log so it shows up
            logging.error(e, exc_info=True)

        # still raise exception
        raise
