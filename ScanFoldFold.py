#!/usr/bin/python3
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

Contact: Ryan Andrews - randrews@iastate.edu

Usage:
$ python3.6 ScanFold-Fold2.0.py fasta_filename [options]

"""


from ScanFoldFunctions import *


import math
import itertools
import operator
from collections import Counter, defaultdict
import sys
import re
import numpy as np
import os
#import RNAstructure
import time
import argparse
from itertools import repeat
from functools import partial
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import shutil

import logging

def main(args):
    start_time = time.time()

    filename = args.filename
    filter = int(args.f)
    competition = int(args.c)
    name = args.id
    global_refold = args.global_refold
    folder_name = args.folder_name
    temperature = int(args.t)
    algo = "rnafold"
    type = "mono"
    cwd = os.getcwd()
    es_path = args.es_path
    igv_path = args.igv_path
    inforna_path = args.inforna_path

    ###create a folder for extracted_structures:
    # extract_directory = "extracted_structures"
    parent_directory = str(os.getcwd())
    extract_path = os.path.join(parent_directory, es_path)
    os.mkdir(extract_path)

    ###create a folder for igv files:
    # igv_directory = "igv_files"
    parent_directory = str(os.getcwd())
    #print(parent_directory)
    igv_path = os.path.join(parent_directory, igv_path)
    os.mkdir(igv_path)

    ###create a folder for inforna_structures:
    # inforna_directory = "inforna_structures"
    parent_directory = str(os.getcwd())
    inforna_path = os.path.join(parent_directory, inforna_path)
    os.mkdir(inforna_path)



    ### This sets strand to forward (but can be set up as an argument if needed)
    strand = 1

    with open(filename, 'r') as f:
        #Read all lines from ScanFold-Scan file (including header)
        try:
            filename_split = filename.split(".tsv")
            outname = filename_split[0]
        except:
            logging.info("Input name should have .tsv extension")
            outname = filename

        # handle case of input file being an absolute path
        outname = os.path.basename(outname)

        if folder_name != None:
            try:
                logging.info("Putting output in folder named:"+folder_name)
                full_output_path = os.path.join(cwd, folder_name)
                os.mkdir(full_output_path)
                os.chdir(full_output_path)
                folder_name = str(folder_name)
            except:
                logging.info(f"WARNING: FOLDER NAME [{folder_name}] IS NOT UNIQUE")
                raise FileExistsError("Folder name exists. Manually input unique folder name via --folder_name [FOLDER NAME]")
        else:
            try:
                folder_name = outname+".Fold_output"
                logging.info("Putting output in folder named:"+outname)
                full_output_path = os.path.join(cwd, outname+"_fold")
                os.mkdir(full_output_path)
                os.chdir(full_output_path)
                folder_name = str(outname+"_fold")
            except:
                logging.info(f"WARNING: FOLDER NAME [{folder_name}] IS NOT UNIQUE")
                raise FileExistsError("Folder name exists. Manually input unique folder name via --folder_name [FOLDER NAME]")

        lines = f.readlines()[1:] # skip header
        steplines1 = int(lines[0].split("\t")[0])
        steplines2 = int(lines[1].split("\t")[0])
        step_size = steplines2  - steplines1

        #Initialize bp dictionary and z-score lists
        log_total = open(outname+".ScanFold.log", 'w')
        log_win = open(outname+".ScanFold.FinalPartners.txt", 'w')

        ### Reset file names ###
        dbn_file_path = outname+".AllDBN_refold.txt"
        dbn_file_path1 = outname+".NoFilter"
        dbn_file_path2 = outname+".Zavg_-1"
        dbn_file_path3 = outname+".Zavg_-2"
        dbn_file_path4 = outname+".AllDBN.txt"
        structure_extract_file = outname+".structure_extracts.txt"

        z_score_list = []
        bp_dict = {}

        #Generate nucleotide dictionary to assign each nucleotide in sequence a key
        nuc_dict = NucleotideDictionary(lines)
        full_fasta_sequence = nuc_dict_to_seq(nuc_dict)
        logging.info("Sequence length: "+str(len(nuc_dict))+"nt")
        # if len(nuc_dict) > 1000:
        #     raise SystemExit('Input sequence is longer than 1000 nt; in order to scan longer sequences consider using the stand alone programs (avaiable here: https://github.com/moss-lab/ScanFold)')

        #Determine start and end coordinate values
        start_coordinate = str(list(nuc_dict.keys())[0])
        #logging.info(start_coordinate)
        end_coordinate = str(list(nuc_dict.keys())[-1])
        #logging.info(end_coordinate)

        #initiate value lists
        MFE_list = []
        zscore_list = []
        pvalue_list = []
        ED_list = []

        #Iterate through input file, read each rows metrics, sequence, etc.
        logging.info("Reading sequence and structures...")
        for row in lines:

            #Ignore blank lines
            if not row.strip():
                continue

            #Main loop to find all i-j pairs per i-nucleotide
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

                strand = 1

                sequence = list(sequence_raw)
                structure = list(structure_raw)

                #Define window coordinates as string
                #window = str(str(icoordinate)+"-"+str(jcoordinate))

                #Determine length of window
                length = len(sequence)

                #Append window z-score to list (to calculate overall z-score)
                z_score_list.append(zscore)

                #Iterate through dot bracket structure to determine locations of bps
                i = 0
                while i < length:
                    #Unpaired nucleotide
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
                    #Paired nucleotide
                    else:
                        i += 1

                #Inititate base pair tabulation variables
                bond_order = []
                bond_count = 0

                #Iterate through sequence to assign nucleotides to structure type
                m = 0
                while  m < length:
                    if structure[m] == '(':
                        bond_count += 1
                        bond_order.append(bond_count)
                        m += 1
                    elif structure[m] == ')':
                        bond_order.append(bond_count)
                        bond_count -= 1
                        m += 1
                    elif structure[m] == '.':
                        bond_order.append(0)
                        m += 1
                    else:
                        logging.info("Error")

                #Initiate base_pair list
                base_pairs = []

                #Create empty variable named test
                test = ""

                #Iterate through bond order
                j = 0
                while j < len(bond_order):
                    if bond_order[j] != 0:
                        test = bond_order[j]
                        base_pairs.append(j+1)
                        bond_order[j] = 0
                        j += 1
                        k = 0
                        while k < len(bond_order):
                            if bond_order[k] == test:
                                base_pairs.append(k+1)
                                bond_order[k] = 0
                                k += 1
                            else:
                                k += 1
                    else:
                        j += 1

                #Iterate through "base_pairs" "to define bps
                l = 0
                while l < len(base_pairs):
                    lbp = base_pairs[l]
                    rbp = base_pairs[l+1]

                    lb = str(sequence[int(lbp)-1])
                    rb = str(sequence[int(rbp)-1])

                    lbp_coord = int(int(lbp)+int(icoordinate)-1)
                    rbp_coord = int(int(rbp)+int(icoordinate)-1)
                    x = NucPair(lb, lbp_coord, rb, rbp_coord, zscore, mfe, ed)
                    z = NucPair(rb, rbp_coord, lb, lbp_coord, zscore, mfe, ed)

                    #Try to append i-j pair to i-nuc for left i-nuc
                    try:
                        y = bp_dict[lbp_coord]
                        y.append(x)
                    #If i-nuc not defined, define it
                    except:
                        bp_dict[lbp_coord] = []
                        y = bp_dict[lbp_coord]
                        y.append(x)

                    #Try to append i-j pair to i-nuc for right i-nuc
                    try:
                        w = bp_dict[rbp_coord]
                        w.append(z)
                    #If i-nuc not defined, define it
                    except:
                        bp_dict[rbp_coord] = []
                        w = bp_dict[rbp_coord]
                        w.append(z)
                    l += 2
                MFE_list.append(mfe)
                zscore_list.append(zscore)
                pvalue_list.append(pvalue)
                ED_list.append(ed)
            #Define OVERALL values of metrics
            meanz = float(np.mean(z_score_list))
            sdz = float(np.std(z_score_list))
            minz = min(z_score_list)
            stdz = float(np.std(z_score_list))

            one_sig_below = float(meanz-stdz)
            two_sig_below = float( meanz - ( 2 * stdz) )

    #Initiate global dictionaries to store best base pairs
    best_bps = {}
    best_sum_bps = {}
    best_sum_bps_means = {}
    best_total_window_mean_bps = {}

    #Iterate through initial i-nuc dictionary to determine best base pairs (round 1)
    elapsed_time = round((time.time() - start_time), 2)
    logging.info("Elapsed time: "+str(elapsed_time)+"s")
    logging.info("Determining best base pairs...")
    for k, v in sorted(bp_dict.items()):
        #Initiate local dictionaries to store metrics per nucleotide
        zscore_dict = {}
        pair_dict = {}
        mfe_dict = {}
        ed_dict = {}

        #Iterate through all i-j pairs per i-nucleotide to store metrics for each
        for pair in v:
            #Create a key  which contains nucleotide and coordinate info
            partner_key = str(pair.jnucleotide)+"-"+str(pair.jcoordinate)
            #logging.info(partner_key)

            #Create a variable which contains all i-j pair info
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)

            #Try to append the value of each metric to metric lists per i-nuc
            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)

                m = mfe_dict[partner_key]
                m.append(pair.mfe)

                e = ed_dict[partner_key]
                e.append(pair.ed)

                z = pair_dict[partner_key]
                z.append(x)

            #If pair not defined, define it
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


        #Caclulate and store sum of z-score per i-j pair
        sum_z = {}
        sum_z_lengths = {}
        for k1, v1 in zscore_dict.items():
            sum_z[k1] = sum(v1)
            test = sum_z[k1] = sum(v1)
            sum_z_lengths[k1] = len(sum_z)

        #Caclulate and store mean of z-score per i-j pair
        mean_z = {}
        for k1, v1 in zscore_dict.items():
            mean_z[k1] = np.mean(v1)
            test = mean_z[k1] = np.mean(v1)

        #Caclulate and store mean MFE per i-j pair
        mean_mfe = {}
        for k1, v1 in mfe_dict.items():
            mean_mfe[k1] = np.mean(v1)

        #Caclulate and store mean ED per i-j pair
        mean_ed = {}
        for k1, v1 in ed_dict.items():
            mean_ed[k1] = np.mean(v1)

        #Caclulate and store total window counts per i-j pair
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

        #Print first line of log file tables (first half of log file)
        k_nuc = str((nuc_dict[k].nucleotide))
        #logging.info(k_nuc)
        log_total.write("\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED"
            +"\tSumZ\tSumZ/#TotalWindows\tBPs= "+str(num_bp)+"\n")
        log_total.write("nt-"+str(k)+"\t-\t"+str(k_nuc)+"\t"+str(total_windows)
            +"\t-\t-\t-\t-\t-"+"\n")


        #Print remainder of log file tables (first half of log file)
        total_window_mean_z = {}
        for k1, v1 in zscore_dict.items():
            bp_window = str(len(v1))
            key_data = re.split("-", str(k1))
            key_nuc = str(key_data[0])
            key_i = str(key_data[1])
            total_window_mean_z[k1] = (sum(v1))/total_windows
            z_sum = str(round(sum(v1), 2))
            z_avg = str(round(np.mean(v1), 2))
            test = str(round(total_window_mean_z[k1], 2))
            k1_mean_mfe = str(round(mean_mfe[k1], 2))
            k1_mean_ed = str(round(mean_ed[k1], 2))
            if int(k) == int(key_i):
                #logging.info("iNuc is "+str(key_i))
                log_total.write(str(k)+"\tNoBP\t"+key_nuc+"\t"+bp_window+"\t"
                    +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                    +z_sum+"\t"+str(test)+"\n")
            else:
                #logging.info("j is "+str(k))
                log_total.write(str(k)+"\t"+key_i+"\t"+key_nuc+"\t"+bp_window+"\t"
                    +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                    +z_sum+"\t"+str(test)+"\n")

        #Define best_bp_key based on coverage-normalized z-score
        best_bp_key = min(total_window_mean_z, key = total_window_mean_z.get)

        #Access best i-j NucPairs for each metric using best_bp_key
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

    ######## Detect competing partners, and select final i-j pairs #################
    final_partners = {}
    elapsed_time = round((time.time() - start_time), 2)
    logging.info("Elapsed time: "+str(elapsed_time)+"s")

    #print header for fianl partener log file (log_win)
    log_win.write("\ni\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED"
        + "\t*Indicates most favorable bp has competition; bp(j) has more favorable partner or is "
        + "more likely to be unpaired"+"\n")

    #Iterate through round 1 i-j pairs
    if competition == 1:
        #logging.info(start_coordinate, end_coordinate)
        logging.info("Detecting competing pairs...")
        j_coord_list = []
        # for k, v in sorted(best_bps.items()):
        #     logging.info(jcoordinate)
        #     j_coord_list.append(int(v.jcoordinate))

        for k, v in sorted(best_bps.items()):
            #logging.info(k, v.icoordinate, v.jcoordinate)
            test_k = int(k)
            #logging.info(sum(test_k == int(v.jcoordinate) for v in best_bps.values()))
            if sum(test_k == int(v.jcoordinate) for v in best_bps.values()) >= 0:
                #print(start_coordinate, end_coordinate)
            ### Scan the entire dictionary:
                # keys = range(int(start_coordinate), int(end_coordinate))

            ### Scan two window's length flanking nucleotide:
                #print(len(best_bps))
                # print(length*4)
                length = len(nuc_dict)
                #print(length)
                if (len(best_bps) < length*4):
                    #print("Scanning full dictionary")
                    #Length of input less than length of flanks
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

                    #logging.info("Found competing pair for "+str(k))
                    #elapsed_time = round((time.time() - start_time), 2)
                    #logging.info(elapsed_time)
                    #logging.info("Detecting competing pairs for nuc ", k)
                    #For each i and j in i-j pair, detect competing pairs and append to dict
                    comp_pairs_i = competing_pairs(subdict, v.icoordinate)
                    #logging.info("i-pairs="+str(len(comp_pairs_i)))
                    comp_pairs_j = competing_pairs(subdict, v.jcoordinate)
                    #logging.info("j-pairs="+str(len(comp_pairs_j)))
                    total_pairs = []
                    #Put pairs competing with i from i-j pair into total pair dict for i-nuc
                    for key, pair in comp_pairs_i.items():
                        #logging.info("checking competing pairs for i")
                        # if k == 216:
                        #     logging.info(k, pair.icoordinate, pair.jcoordinate, pair.zscore)
                        total_pairs.append(competing_pairs(subdict,
                                                        pair.jcoordinate))
                        total_pairs.append(competing_pairs(subdict,
                                                        pair.icoordinate))
                    #Put pairs competing with j from i-j pair into total pair dict for i-nuc
                    for key, pair in comp_pairs_j.items():
                        #logging.info("checking competing pairs for j")
                        # if k == 216:
                        #     logging.info(k, pair.icoordinate, pair.jcoordinate, pair.zscore)
                        total_pairs.append(competing_pairs(subdict,
                                                        pair.jcoordinate))
                        total_pairs.append(competing_pairs(subdict,
                                                        pair.icoordinate))
                    # logging.info(str(k)+"nt Total comp pairs="+str(len(total_pairs)))

                    #Merge all dictionaries
                    merged_dict = {}
                    i = 0
                    for d in total_pairs:
                        #logging.info("merging competing dictionaries "+str(i))
                        for k1, v1 in d.items():
                            # if k == 216:
                            #     logging.info(k, k1, v1.icoordinate, v1.jcoordinate, v1.zscore)
                            merged_dict[i] = v1
                            i += 1

                    # #logging.info("MergedDict length for "+str(k)+"="+str(len(merged_dict)))
                    # #initiate best_basepair fucntion, return best_bp based on sum
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
                        bp = NucPair(v.inucleotide, v.icoordinate,
                                                    v.inucleotide, v.icoordinate,
                                                    v.zscore,
                                                    v.mfe,
                                                    v.ed)

                    if (int(k) != bp.icoordinate) and (int(k) != int(bp.jcoordinate)):
                        # logging.info("1 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                        #if there was a competing i-j pair print it to log file instead:
                        log_win.write("nt-"+str(k)+"*:\t"+str(v.icoordinate)+"\t"+v.jcoordinate+"\t"
                            +str(round(v.mfe, 2))
                            +"\t"+str(round(v.zscore, 2))
                            +"\t"+str(round(v.ed, 2))+"\n")
                        final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                                    v.inucleotide, v.icoordinate,
                                                    best_bps[bp.icoordinate].zscore,
                                                    bp.mfe,
                                                    bp.ed)
                    #
                    # elif ((int(v.icoordinate) == int(v.jcoordinate)) and (int(bp.icoordinate) != int(bp.jcoordinate))):
                    #     #Check for instance where competing base pair
                    #     logging.info("!!!!!!!2 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                    #     log_win.write("nt-"+str(k)+"*:\t"+str(bp.icoordinate)+"\t"+bp.jcoordinate+"\t"
                    #           +str(round(bp.mfe, 2))
                    #           +"\t"+str(round(bp.zscore, 2))
                    #           +"\t"+str(round(bp.ed, 2))+"\n")
                    #
                    #     final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                    #                                 bp.jnucleotide, bp.jcoordinate,
                    #                                 best_bps[bp.icoordinate].zscore,
                    #                                 best_bps[bp.icoordinate].mfe,
                    #                                 best_bps[bp.icoordinate].ed)
                    #
                    #
                    else:
                        #logging.info("3 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                        log_win.write("nt-"+str(k)+":\t"+str(bp.icoordinate)+"\t"+str(bp.jcoordinate)+"\t"
                            + str(round(best_bps[k].mfe, 2))+"\t"
                            + str(round(best_bps[k].zscore, 2))
                            + "\t"+str(round(best_bps[k].ed, 2))+"\n")
                        final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                                                    bp.jnucleotide, bp.jcoordinate,
                                                    best_bps[bp.icoordinate].zscore,
                                                    best_bps[bp.icoordinate].mfe,
                                                    best_bps[bp.icoordinate].ed)

                else:
                    continue
            else:
                final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                            v.jnucleotide, v.jcoordinate,
                                            best_bps[k].zscore,
                                            best_bps[k].mfe,
                                            best_bps[k].ed)
                #logging.info("No competing pair found for ", k)
                continue

    #Write CT files
    if competition == 1:
        logging.info("Trying to write CT files with -c option")
        elapsed_time = str(round((time.time() - start_time), 2))+"s"
        logging.info("Elapsed time: "+elapsed_time)
        logging.info("Writing CT files")

        write_ct(final_partners, str(dbn_file_path1)+".ct", float(10), strand, name, start_coordinate)
        write_ct(final_partners, str(dbn_file_path2)+".ct", float(-1), strand, name, start_coordinate)
        write_ct(final_partners, str(dbn_file_path3)+".ct", float(-2), strand, name, start_coordinate)
        if filter != -2:
            write_ct(final_partners, str(user_filter_dbn)+".ct", filter, strand, name, start_coordinate)
            makedbn(user_filter_dbn, "Zavg"+str(filter))
        makedbn(dbn_file_path1, "NoFilter")
        makedbn(dbn_file_path2, "Zavg_-1")
        makedbn(dbn_file_path3, "Zavg_-2")

        write_bp(final_partners, igv_path+"/"+outname+".bp", start_coordinate, name, minz)
        write_wig_dict(final_partners, igv_path+"/"+outname+".zavgs.wig", name, step_size, str("zscore"))
        write_wig_dict(final_partners, igv_path+"/"+outname+".mfe_avgs.wig", name, step_size, str("mfe"))
        write_wig_dict(final_partners, igv_path+"/"+outname+".ed_avgs.wig", name, step_size, str("ed"))
        write_bp(best_bps, igv_path+"/"+outname+".ALL.bp", start_coordinate, name, minz)


        write_bp(final_partners, outname+".bp", start_coordinate, name, minz)
        if args.bp_track is not None:
            write_bp(final_partners, args.bp_track, start_coordinate, name, minz)

        write_wig_dict(final_partners, outname+".zavgs.wig", name, step_size, str("zscore"))
        if args.final_partners_wig is not None:
            write_wig_dict(final_partners, args.final_partners_wig, name, step_size, str("zscore"))

        write_wig_dict(final_partners, outname+".mfe_avgs.wig", name, step_size, str("mfe"))
        if args.mfe_wig_file_path is not None:
            write_wig_dict(final_partners, args.mfe_wig_file_path, name, step_size, str("mfe"))

        write_wig_dict(final_partners, outname+".ed_avgs.wig", name, step_size, str("ed"))
        if args.ed_wig_file_path is not None:
            write_wig_dict(final_partners, args.ed_wig_file_path, name, step_size, str("ed"))

        write_bp(best_bps, outname+".ALL.bp", start_coordinate, name, minz)


    if competition == 0:
        if args.bp_track is not None:
            write_bp(best_bps, args.bp_track, start_coordinate, name, minz)

    write_fasta(nuc_dict, outname+".fa", name)
    if args.fasta_file_path is not None:
        write_fasta(nuc_dict, args.fasta_file_path, name)

    write_fai(nuc_dict, outname+".fai", name)
    if args.fasta_index is not None:
        write_fai(nuc_dict, args.fasta_index, name)

    logging.info("ScanFold-Fold analysis complete! Refresh page to ensure proper loading of IGV")
    merge_files(str(dbn_file_path4), str(dbn_file_path1+".dbn"), str(dbn_file_path2+".dbn"), str(dbn_file_path3+".dbn"))

    if global_refold == True:
        logging.info("The global_refold section of code is a work in progress.")
        #create a separate DBN file
        dbn_log_file = open(str(dbn_file_path), "w+")

        #try:
        #fold the full fasta input as a fold compound (full_fc) using model params (md)
        logging.info("Refolding full sequence using ScanFold results as constraints...")
        elapsed_time = round((time.time() - start_time), 2)
        logging.info("Elapsed time: "+str(elapsed_time)+"s")
        md = RNA.md()
        md.temperature = int(temperature)

        #refold from -1 constraints
        fc = RNA.fold_compound(str(full_fasta_sequence), md)
        dbn_file_filter1 = open(str(dbn_file_path2+".dbn"), "r")
        lines = dbn_file_filter1.readlines()
        filter1constraints = str(lines[2])
        refolded1filter = fc.hc_add_from_db(filter1constraints)
        (refolded_filter1_structure, refolded_filter1_MFE) = fc.mfe()

        #refold from -2 constraints
        fc = RNA.fold_compound(str(full_fasta_sequence), md)
        dbn_file_filter2 = open(str(dbn_file_path3+".dbn"), "r")
        lines = dbn_file_filter2.readlines()
        filter2constraints = str(lines[2])
        refolded2filter = fc.hc_add_from_db(filter2constraints)
        (refolded_filter2_structure, refolded_filter2_MFE) = fc.mfe()

        #extract the structure
        full_fc = RNA.fold_compound(str(full_fasta_sequence), md)
        (full_structure, full_MFE) = full_fc.mfe()

        dbn_log_file.write(">"+str(name)+"\tGlobal Full MFE="+str(full_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(full_structure)+"\n")
        dbn_log_file.write(">"+str(name)+"\tRefolded with -1 constraints MFE="+str(refolded_filter1_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(refolded_filter1_structure)+"\n")
        dbn_log_file.write(">"+str(name)+"\tRefolded with -2 constraints MFE="+str(refolded_filter2_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(refolded_filter2_structure)+"\n")
        dbn_log_file.close()
        merge_files(str(dbn_file_path4), str(dbn_file_path1+".dbn"), str(dbn_file_path2+".dbn"), str(dbn_file_path3+".dbn"))
        dbn_log_file.close()
        # except:
        #     logging.info("Automatic refold for "+cur_record.name+" failed. Run manually")
        #     os.chdir(str(original_directory))
        #     del final_partners
        #     continue

    if global_refold == False:
        merge_files(str(dbn_file_path4), str(dbn_file_path1+".dbn"), str(dbn_file_path2+".dbn"), str(dbn_file_path3+".dbn"))

    #############
    #Begin the structure extract process
    #Set flanking nucleotides to be folded
    flanking = 0

    #Set number of randomizations and shuffle type ("mono" or "di")
    sub_randomizations = 100
    sub_shuffle = "mono"


    #Inititate base pair tabulation variables
    bond_order = []
    bond_count = 0

    #Read the structure of -2 filter2constraints
    if global_refold == False:
        # Refold from -1 constraints
        if args.extract == 1:
            dbn_file_filter1 = open(dbn_file_path2+".dbn", "r")
            lines = dbn_file_filter1.readlines()
            full_fasta_sequence = str(lines[1])
            filter2constraints = str(lines[2])
        if args.extract == 2:
            dbn_file_filter2 = open(dbn_file_path3+".dbn", "r")
            lines = dbn_file_filter2.readlines()
            full_fasta_sequence = str(lines[1])
            filter2constraints = str(lines[2])

    structure_raw = filter2constraints
    sequence = list(str(full_fasta_sequence))
    structure = list(structure_raw)
    length = len(sequence)
    length_st = len(structure)

    if length != length_st:
        raise("Length of sequence and structure do not match")

    #Iterate through sequence to assign nucleotides to structure type
    m = 0
    nuc_dict_refold = {}
    while  m < len(structure)-1:
        #logging.info(m)
        if structure[m] == '(':
            #logging.info(m, structure[m])
            bond_count += 1
            bond_order.append(bond_count)
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1

        elif structure[m] == ')':
            # logging.info(m, structure[m])
            bond_order.append(bond_count)
            bond_count -= 1
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1

        elif str(structure[m]) == ( '.' ):
            # logging.info(m, structure[m])
            bond_order.append(0)
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '<' ):
            # logging.info(m, structure[m])
            bond_order.append(0)
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '>' ):
            # logging.info(m, structure[m])
            bond_order.append(0)
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '{' ):
            # logging.info(m, structure[m])
            bond_order.append(0)
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '}' ):
            # logging.info(m, structure[m])
            bond_order.append(0)
            nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        else:
            # logging.info("Error", bond_count, (m+1), sequence[m], structure[m])
            m += 1
            continue
            # logging.info("no")

        """ Repeat the process looking for non-nested "<..>" pairs """
        #Inititate base pair tabulation variables
        bond_order_pk = []
        bond_count_pk = 0
        nuc_dict_pk = {}
        #Iterate through sequence to assign nucleotides to structure type
        m = 0
        while  m < length-1:
            #print(m)
            if structure[m] == '<':
                #print(m, structure[m])
                bond_count_pk += 1
                bond_order_pk.append(bond_count_pk)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1

            elif structure[m] == '>':
            #    print(m, structure[m])
                bond_order_pk.append(bond_count_pk)
                bond_count_pk -= 1
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1

            elif str(structure[m]) == ( '.' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '(' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( ')' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '{' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '}' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            else:
                #print("Error", bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
                continue


    #Initiate base_pair list
    base_pairs = []

    #Create empty variable named test
    test = ""

    #Iterate through bond order
    j = 0
    structure_count = 0
    structure_end = []
    structure_start = []
    while j < len(structure)-1:
        #logging.info(nuc_dict[j].structure)
        try:
            if (nuc_dict_refold[j].bond_order == 1) and (nuc_dict_refold[j].structure == '('):
                structure_count += 1
                #print(nuc_dict[j].structure, j)
                structure_start.append(NucStructure(structure_count, nuc_dict_refold[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                j += 1

            elif (nuc_dict_refold[j].bond_order == 0) and (nuc_dict_refold[j].structure == ')'):
                structure_count += 1
                #print(nuc_dict[j].structure, j)
                structure_end.append(NucStructure(structure_count, nuc_dict_refold[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                j += 1
            else:
                j += 1
        except:
            j += 1
            continue

    """ Repeat for non-nested """
    j = 0
    structure_count_pk = 0
    structure_end_pk = []
    structure_start_pk = []

    while j < length:
        #print(f"{j} non-nested")
        try:
            if (nuc_dict_pk[j].bond_order == 1) and (nuc_dict_pk[j].structure == '<'):
                structure_count_pk += 1
                #print(nuc_dict_pk[j].structure)
                structure_start_pk.append(NucStructure(structure_count_pk, nuc_dict_pk[j].coordinate, nuc_dict_pk[j].nucleotide, nuc_dict_pk[j].structure))
                j += 1

            elif (nuc_dict_pk[j].bond_order == 0) and (nuc_dict_pk[j].structure == '>'):
                structure_count_pk += 1
                #print(nuc_dict_pk[j].structure)
                structure_end_pk.append(NucStructure(structure_count, nuc_dict_pk[j].coordinate, nuc_dict_pk[j].nucleotide, nuc_dict_pk[j].structure))
                j += 1
            else:
                j += 1
        except:
            #print("Here")
            j += 1
            continue

    """Add to extracted structure list"""
    l = 0
    extracted_structure_list = []
    #logging.info(structure_start)
    #logging.info(structure_end)
    while l < int(len(structure_start)):
        #logging.info(len(structure_start), len(structure_end))
        offset = flanking = 0
        s = structure_start_coordinate =  int((structure_start[l].coordinate)-offset)
        e = structure_end_coordinate = int((structure_end[l].coordinate)+offset-1)

        seq = ""
        se_fold = ""
        for k, v in nuc_dict_refold.items():
            if s <= k <= e:
                seq += str(v.nucleotide)
                se_fold += str(v.structure)

        extracted_structure_list.append(ExtractedStructure(l, seq, se_fold, s, e))
        l += 1

    """ Repeat for non-nested """
    l = 0
    while l < int(len(structure_start_pk)):
        offset = flanking
        s = structure_start_coordinate =  int((structure_start_pk[l].coordinate)-offset-1)
        e = structure_end_coordinate = int((structure_end_pk[l].coordinate)+offset-1)

        seq = ""
        fold = ""
        for k, v in nuc_dict_pk.items():
            if s <= k <= e:
                seq += str(v.nucleotide)
                fold += str(v.structure)

        extracted_structure_list.append(ExtractedStructure(l, seq, fold, s, e))

        l += 1

    zscore_total = []
    numerical_z = []
    pvalue_total = []
    numerical_p = []
    MFE_total = []
    ED_total = []

    se = open(structure_extract_file, "w")
    #se.write("ScanFold predicted structures which contain at least one base pair with Zavg < -1 have been extracted from "+str(name)+" results (sequence length "+str(length)+"nt) and have been refolded using RNAfold to determine their individual MFE, structure, z-score (using 100X randomizations), and ensemble diversity score.\n")
    motif_num = 1
    for es in extracted_structure_list[:]:
        frag = es.sequence
        fc = RNA.fold_compound(str(frag)) #creates "Fold Compound" object
        fc.hc_add_from_db(str(es.structure))
        fc.pf() # performs partition function calculations
        frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
        (MFE_structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
        MFE = round(MFE, 2)
        MFE_total.append(MFE)
        (centroid, distance) = fc.centroid() # calculate and define variables for centroid
        ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
        ED_total.append(ED)
        seqlist = [] # creates the list we will be filling with sequence fragments
        seqlist.append(frag) # adds the native fragment to list
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
        ps_title = f"motif_{motif_num} coordinates {es.i} - {es.j}"
        es_path = f"{extract_path}/{name}_motif_{motif_num}.dbn"
        with open(es_path, 'w') as es_dbn:
            es_dbn.write(f">{name}_motif_{motif_num}_coordinates:{es.i}-{es.j}_zscore={zscore}\n{frag}\n{MFE_structure}")
        dbn2ct(es_path)
        os.rename(f"{extract_path}/{name}_motif_{motif_num}.ct", f"{inforna_path}/{name}_motif_{motif_num}.ct")

        # Create postcript files
        RNA.PS_rna_plot_a(frag, MFE_structure, extract_path+"/motif_"+str(motif_num)+".ps", '', '')

        # Set extracted structures up as GFF format
        gff_attributes = f'motif_{motif_num};sequence={es.sequence};structure={str(es.structure)};refoldedMFE={str(MFE_structure)};MFE(kcal/mol)={str(MFE)};z-score={str(zscore)};ED={str(ED)}'
        se.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str("."), str("RNA_sequence_secondary_structure"), str(int((es.i)+1)), str(int((es.j)+1)), str("."), str("."), str("."), gff_attributes))
        motif_num += 1

    se.close()
    # except:
    #     logging.info("Structure Extract failed for "+folder_name+", must extract manually.")
    #     continue
    #logging.info("TEST")
    shutil.make_archive(inforna_path, 'zip', inforna_path)
    elapsed_time = round((time.time() - start_time), 2)
    logging.info("Total runtime: "+str(elapsed_time)+"s")
    logging.info("ScanFold-Fold analysis complete! Output found in folder named: "+folder_name)

    if args.webserver:
        make_tar(args.webserver, full_output_path)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('filename',  type=str,
                        help='input filename')
    parser.add_argument('-f', type=int, default=-2,
                        help='filter value')
    parser.add_argument('-c', type=int, default=1,
                        help='Competition (1 for disallow competition, 0 for allow; 1 by default)')
    parser.add_argument('--id', type=str, default = "UserInput",
                        help='Name or ID of sequence being analyzed. Default "UserInput"')
    parser.add_argument('--global_refold', action='store_true', default=False,
                        help='Global refold option. Refold full sequence using Zavg <-1 and <-2 base pairs')
    parser.add_argument('--folder_name',  type=str,
                        help='Name of output folder (defaults to header name or date/time)', default=None)
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')
    parser.add_argument('--extract', type=int, default='2',
                        help='Extract structures from minus 1 or minus 2 dbn file (2 or 1); Default = 2')
    parser.add_argument('--es_path', type=str, default = "extracted_structures",
                        help='')
    parser.add_argument('--igv_path', type=str, default = "igv_files",
                        help='')
    parser.add_argument('--inforna_path', type=str, default = "inforna_structures",
                        help='')


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

    logging.basicConfig(stream=args.logfile, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=loglevel)

    try:
        main(args)
    except Exception as e:

        if args.webserver:
            # log so it shows up
            logging.error(e, exc_info=True)

        # still raise exception
        raise
