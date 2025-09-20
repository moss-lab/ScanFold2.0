#!/usr/bin/python3
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

ScanFold2.0

Contact: Ryan Andrews - randrews@iastate.edu

Usage:
$ python3.6 ScanFold.py fasta_filename [options]

"""
import time
start_time = time.time()

from ScanFoldFunctions import *
import argparse
from Bio import SeqIO
import pandas as pd
import tensorflow as tf
from datetime import datetime
import os
import sys
import json

### Arguments
parser = argparse.ArgumentParser()
parser.add_argument('filename',  type=str,
                    help='input filename')
parser.add_argument('-t', type=int, default=37,
                    help='Folding temperature in celsius; default = 37C')
parser.add_argument('-s', type=int, default=1,
                    help='Step size; default = 1')
parser.add_argument('-w', type=int, default=120,
                    help='Window size; default = 120')
parser.add_argument('--folder', type=str, default=None,
        help='Output folder name (default: ScanFold_YYYYMMDD_HHMMSS)')
parser.add_argument('--shuffle', type=str, default='mono', choices=['mono', 'di'],
        help='Shuffling type: mono (mononucleotide) or di (dinucleotide); default = mono')
#parser.add_argument('-a', type=str, default='rnafold',
#        help='Folding algorithm to use; default = rnafold')

### Parse arguments and convert to variables
args = parser.parse_args()
myfasta = args.filename
temperature = int(args.t)
step_size = int(args.s)
window_size = int(args.w)
output_folder = args.folder
shuffle_type = args.shuffle
#algo = str(args.a)

# Get absolute path of input file BEFORE changing directories
myfasta = os.path.abspath(myfasta)
basename = os.path.basename(myfasta).split('.')[0]

### Create output folder
if output_folder is None:
    # Generate timestamp-based folder name
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_folder = f"ScanFold_{timestamp}"

# Create the output folder if it doesn't exist (use absolute path)
output_folder = os.path.abspath(output_folder)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created output folder: {output_folder}")

# Store original directory for model loading
original_dir = os.getcwd()

# Create log file with absolute path
log_file = os.path.join(output_folder, "ScanFold_run.log")
with open(log_file, 'w') as log:
    log.write("="*60 + "\n")
    log.write("ScanFold2.0 Run Log\n")
    log.write("="*60 + "\n")
    log.write(f"Date/Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"Input file: {myfasta}\n")
    log.write(f"Output folder: {output_folder}\n")
    log.write(f"\nParameters:\n")
    log.write(f"  Temperature: {temperature}°C\n")
    log.write(f"  Window size: {window_size} nt\n")
    log.write(f"  Step size: {step_size} nt\n")
    log.write(f"  Shuffle type: {shuffle_type}\n")
    log.write(f"  Python version: {sys.version.split()[0]}\n")
    log.write(f"  TensorFlow version: {tf.__version__}\n")
    log.write("\n" + "="*60 + "\n\n")

# Function to append to log file
def log_message(message):
    with open(log_file, 'a') as log:
        timestamp = datetime.now().strftime('%H:%M:%S')
        log.write(f"[{timestamp}] {message}\n")
        print(message)

# Helper function to get output file path
def get_output_path(filename):
    return os.path.join(output_folder, filename)

log_message(f"Starting ScanFold analysis...")

### Load appropriate models based on shuffle type
# Use local model paths relative to script directory (get before changing dir)
script_dir = os.path.dirname(os.path.abspath(__file__))

if shuffle_type == 'di':
    # Dinucleotide shuffling models
    log_message("Using dinucleotide shuffling models (DiMFE/DiStd)")
    mfe_model = tf.keras.models.load_model(os.path.join(script_dir, 'DiMFE'))
    stddev_model = tf.keras.models.load_model(os.path.join(script_dir, 'DiStd'))
else:
    # Mononucleotide shuffling models (default)
    log_message("Using mononucleotide shuffling models (MeanMFE/StdDev)")
    mfe_model = tf.keras.models.load_model(os.path.join(script_dir, 'MeanMFE'))
    stddev_model = tf.keras.models.load_model(os.path.join(script_dir, 'StdDev'))

### Start main loop
with open(myfasta, 'r') as forward_fasta:
    for cur_record in SeqIO.parse(forward_fasta, "fasta"):
        df = pd.DataFrame(columns = ["Start", "End", "Temperature", "NativeMFE",
            "Z-score", "p-value", "ED", "Sequeunce", "Structure", "centroid"])
        read_name = cur_record.name
        name = read_name
        ### Create output files with absolute paths
        output = str(read_name+"."+os.path.basename(myfasta)+".ScanFold.")
        outname = str(read_name+".win_"+str(window_size)+".stp_"+str(step_size)+".csv")
        log_total = open(get_output_path(outname+".ScanFold.log"), 'w')
        log_win = open(get_output_path(outname+".ScanFold.FinalPartners.txt"), 'w')

        ### Grab sequence and its features
        cur_record_length = len(cur_record.seq)
        seq = cur_record.seq
        seq = seq.transcribe()

        ### Create a nucleotide dictionary
        nuc_dict = {}
        nuc_i = 1
        for nuc in seq:
            x = NucZscore(nuc,(nuc_i))
            nuc_dict[x.coordinate] = x
            nuc_i += 1

        (start_coordinate_list, end_coordinate_list, frag_list,
         length_list, GCpercent_list, CGratio_list, AUratio_list,
         mfe_list, structure_list, centroid_list,
         ed_list, di_freq_list) = get_frag_feature_list(seq, step_size, window_size, 'rnafold', temperature)

        five_feature_predict = pd.DataFrame(columns = ["Length", "GCpercent","CGratio", "AUratio", "MFE"])
        five_feature_predict["Length"] = length_list
        five_feature_predict["GCpercent"] = GCpercent_list
        five_feature_predict["CGratio"] = CGratio_list
        five_feature_predict["AUratio"] = AUratio_list
        five_feature_predict["MFE"] = mfe_list

        four_feature_predict = pd.DataFrame(columns = ["Length", "GCpercent","CGratio", "AUratio"])
        four_feature_predict["Length"] = length_list
        four_feature_predict["GCpercent"] = GCpercent_list
        four_feature_predict["CGratio"] = CGratio_list
        four_feature_predict["AUratio"] = AUratio_list

        # Combine with dinucleotide frequencies to create full feature set
        full_feature_predict = four_feature_predict
        di_freq_df = pd.DataFrame(di_freq_list, columns = ["AA","AU","AG","AC","UA","UU","UG","UC","GA","GU","GG","GC","CA","CU","CG", "CC"])
        full_feature_predict = full_feature_predict.join(di_freq_df)

        meanMFE_result = mfe_model.predict(full_feature_predict)
        stddev_result = stddev_model.predict(full_feature_predict)

        zscore_list = []
        for i in range(0, len(mfe_list)):
            zscore = round(float((mfe_list[i]-meanMFE_result[i][0])/stddev_result[i][0]), 2)
            zscore_list.append(zscore)

        df["Start"] = start_coordinate_list
        df["End"] = end_coordinate_list
        df["Temperature"] = list([temperature] * len(mfe_list))
        df["NativeMFE"] = mfe_list
        df["Z-score"] = zscore_list
        df["p-value"] = list([0] * len(mfe_list))
        df["ED"] = ed_list
        df["Sequence"] = frag_list
        df["Structure"] = structure_list
        df["centroid"] = centroid_list
        df.to_csv(get_output_path(outname), sep="\t")

        minz = df["Z-score"].min()
        log_message(f"Minimum Z-score: {minz}")
        elapsed_time = round((time.time() - start_time), 2)
        log_message(f"ScanFold-Scan complete. Elapsed time: {elapsed_time}s")

        """
        This is the beginning of "ScanFold-Fold".
        """
        bp_dict = {}

        # Vectorized processing of unpaired nucleotides
        # Convert to numpy arrays for faster processing
        import numpy as np

        # Process all rows at once instead of iterating
        sequences = df["Sequence"].values
        structures = df["Structure"].values
        starts = df["Start"].values.astype(int)
        zscores = df["Z-score"].values
        mfes = df["NativeMFE"].values
        eds = df["ED"].values

        # Process each window efficiently
        for idx in range(len(df)):
            sequence = sequences[idx]
            structure = structures[idx]
            start = starts[idx]
            zscore = zscores[idx]
            mfe = mfes[idx]
            ed = eds[idx]

            # Use numpy to find unpaired positions
            unpaired_positions = np.array([i for i, s in enumerate(structure) if s == '.'])

            # Batch create NucPair objects for unpaired nucleotides
            for pos in unpaired_positions:
                nucleotide = sequence[pos]
                coordinate = pos + start
                x = NucPair(nucleotide, coordinate, nucleotide, coordinate,
                           zscore, mfe, ed)
                if coordinate not in bp_dict:
                    bp_dict[coordinate] = []
                bp_dict[coordinate].append(x)

            # Optimized base pair finding using stack-based approach
            base_pairs = []
            stack = []

            # Find base pairs efficiently with a single pass
            for i, char in enumerate(structure):
                if char == '(':
                    stack.append(i + 1)  # Store 1-based position
                elif char == ')':
                    if stack:
                        left_pos = stack.pop()
                        base_pairs.extend([left_pos, i + 1])  # Store both positions

            # Process base pairs efficiently
            for l in range(0, len(base_pairs), 2):
                lbp = base_pairs[l]
                rbp = base_pairs[l+1]

                lb = sequence[lbp-1]
                rb = sequence[rbp-1]

                lbp_coord = lbp + start - 1
                rbp_coord = rbp + start - 1
                x = NucPair(lb, lbp_coord, rb, rbp_coord, zscore, mfe, ed)
                z = NucPair(rb, rbp_coord, lb, lbp_coord, zscore, mfe, ed)

                # Use dict.setdefault for cleaner code
                bp_dict.setdefault(lbp_coord, []).append(x)
                bp_dict.setdefault(rbp_coord, []).append(z)

        best_bps = {}
        best_sum_bps = {}
        best_sum_bps_means = {}
        best_total_window_mean_bps = {}

        #Determine start and end coordinate values
        start_coordinate = str(list(nuc_dict.keys())[0])
        end_coordinate = str(list(nuc_dict.keys())[-1])
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

            #Caclulate and store sum of ED per i-j pair
            sum_ed = {}
            sum_ed_lengths = {}
            for k1, v1 in zscore_dict.items():
                sum_ed[k1] = sum(v1)
                test = sum_ed[k1] = sum(v1)
                sum_ed_lengths[k1] = len(sum_ed)

            #Caclulate and store mean of z-score per i-j pair
            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = statistics.mean(v1)
                test = mean_z[k1] = statistics.mean(v1)

            #Caclulate and store mean MFE per i-j pair
            mean_mfe = {}
            for k1, v1 in mfe_dict.items():
                mean_mfe[k1] = statistics.mean(v1)

            #Caclulate and store mean ED per i-j pair
            mean_ed = {}
            for k1, v1 in ed_dict.items():
                mean_ed[k1] = statistics.mean(v1)

            #Caclulate and store total window counts per i-j pair
            total_windows = 0
            num_bp = 0
            for k1, v1 in zscore_dict.items():
                total_windows = total_windows + len(v1)
                #key_data = re.split("-", str(k1))
                key_i = (k1[2:])
                if int(k) == int(key_i):
                    continue
                if int(k) != int(key_i):
                    num_bp += 1

            #Print first line of log file tables (first half of log file)
            k_nuc = str((nuc_dict[k].nucleotide))
            log_total.write("\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED"
                  +"\tSumZ\tSumZ/#TotalWindows\tBPs= "+str(num_bp)+"\n")
            log_total.write("nt-"+str(k)+"\t-\t"+str(k_nuc)+"\t"+str(total_windows)
                  +"\t-\t-\t-\t-\t-"+"\n")

            #Print remainder of log file tables (first half of log file)
            total_window_mean_z = {}
            for k1, v1 in zscore_dict.items():
                bp_window = str(len(v1))
                key_nuc = k1[:1]
                key_i = int(k1[2:])
                total_window_mean_z[k1] = (sum(v1))/total_windows
                z_sum = str(round(sum(v1), 2))
                z_avg = str(round(statistics.mean(v1), 2))
                test = str(round(total_window_mean_z[k1], 2))
                k1_mean_mfe = str(round(mean_mfe[k1], 2))
                k1_mean_ed = str(round(mean_ed[k1], 2))
                if int(k) == int(key_i):
                    log_total.write(str(k)+"\tNoBP\t"+key_nuc+"\t"+bp_window+"\t"
                          +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                          +z_sum+"\t"+str(test)+"\n")
                else:
                    log_total.write(str(k)+"\t"+str(key_i)+"\t"+str(key_nuc)+"\t"+bp_window+"\t"
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

        #print header for fianl partener log file (log_win)
        log_win.write("i\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED"
            + "\t*Indicates most favorable bp has competition; bp(j) has more favorable partner or is "
            + "more likely to be unpaired"+"\n")

        #Iterate through round 1 i-j pairs
        elapsed_time = round((time.time() - start_time), 2)
        log_message(f"ScanFold-Fold step 1 of 2 completed. Elapsed time: {elapsed_time}s")
        log_message("Detecting competing pairs...")
        j_coord_list = []
        for k, v in sorted(best_bps.items()):
            test_k = int(k)
            if sum(test_k == int(v.jcoordinate) for v in best_bps.values()) >= 0:
                ### Scan two window's length flanking nucleotide:
                length = len(cur_record.seq)
                if (len(best_bps) < length*4):
                    subdict = best_total_window_mean_bps
                elif (
                    (v.icoordinate - length*(2)) >= int(start_coordinate) and
                    (v.icoordinate + (length*2)) <= int(end_coordinate)
                    ):
                    keys = range(int(v.icoordinate-(length*2)), int(v.icoordinate+(length*2)))
                    subdict = {k: best_total_window_mean_bps[k] for k in keys}
                elif (
                    int(v.icoordinate + (length*(2))) <= int(end_coordinate)and
                    (v.icoordinate + (length*2)) <= int(end_coordinate)
                    ):
                    keys = range(int(start_coordinate), int(v.icoordinate+(length*2))+1)
                    subdict = {k: best_total_window_mean_bps[k] for k in keys}
                elif (v.icoordinate + (length*2)) >= int(end_coordinate):
                    if v.icoordinate-(length*2) > 0:
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
                    #For each i and j in i-j pair, detect competing pairs and append to dict
                    comp_pairs_i = competing_pairs(subdict, v.icoordinate)
                    comp_pairs_j = competing_pairs(subdict, v.jcoordinate)
                    total_pairs = []

                    #Put pairs competing with i from i-j pair into total pair dict for i-nuc
                    for key, pair in comp_pairs_i.items():
                        total_pairs.append(competing_pairs(subdict,
                                                           pair.jcoordinate))
                        total_pairs.append(competing_pairs(subdict,
                                                           pair.icoordinate))

                    #Put pairs competing with j from i-j pair into total pair dict for i-nuc
                    for key, pair in comp_pairs_j.items():
                        total_pairs.append(competing_pairs(subdict,
                                                           pair.jcoordinate))
                        total_pairs.append(competing_pairs(subdict,
                                                           pair.icoordinate))

                    #Merge all dictionaries
                    merged_dict = {}
                    i = 0
                    for d in total_pairs:
                        for k1, v1 in d.items():
                            merged_dict[i] = v1
                            i += 1

                    if len(merged_dict) > 0:
                        bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                    else:
                        bp = NucPair(v.inucleotide, v.icoordinate,
                                                    v.inucleotide, v.icoordinate,
                                                    v.zscore,
                                                    v.mfe,
                                                    v.ed)

                    if (int(k) != bp.icoordinate) and (int(k) != int(bp.jcoordinate)):
                        log_win.write("nt-"+str(k)+"*:\t"+str(v.icoordinate)+"\t"+v.jcoordinate+"\t"
                              +str(round(v.mfe, 2))
                              +"\t"+str(round(v.zscore, 2))
                              +"\t"+str(round(v.ed, 2))+"\n")
                        final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                                    v.inucleotide, v.icoordinate,
                                                    best_bps[bp.icoordinate].zscore,
                                                    bp.mfe,
                                                    bp.ed)
                    else:
                        #print("3 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
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
                continue
        elapsed_time = round((time.time() - start_time), 2)
        log_message(f"ScanFold-Fold step 2 of 2 completed. Elapsed time: {elapsed_time}s")
        log_message("Writing output files...")

        # Close analysis log files before writing final outputs
        log_total.close()
        log_win.close()

        strand = 1  # 1 for forward strand, consistent with ScanFoldFold.py
        write_bp(best_bps, get_output_path(basename+cur_record.name+".ALL.bp"), start_coordinate, name, minz)
        write_ct(final_partners, get_output_path(str(basename)+".no_filter.ct"), float(10), strand, name, start_coordinate)
        write_ct(final_partners, get_output_path(str(basename)+".minus_1.ct"), float(-1), strand, name, start_coordinate)
        write_ct(final_partners, get_output_path(str(basename)+".minus_2.ct"), float(-2), strand, name, start_coordinate)
        # makedbn expects base filename without extension and adds .ct/.dbn itself
        # We need to use the full path without extension
        makedbn(os.path.join(output_folder, str(basename)+".no_filter"), "NoFilter")
        makedbn(os.path.join(output_folder, str(basename)+".minus_1"), "Zavg_-1")
        makedbn(os.path.join(output_folder, str(basename)+".minus_2"), "Zavg_-2")
        write_bp(final_partners, get_output_path(outname+".bp"), start_coordinate, name, minz)
        write_wig_dict(final_partners, get_output_path(outname+".Zavg.wig"), name, step_size, "zscore")
        write_wig(mfe_list, step_size, cur_record.name, get_output_path(outname+".scan-MFE.wig"))
        write_wig(zscore_list, step_size, cur_record.name, get_output_path(outname+".scan-zscores.wig"))
        write_wig(ed_list, step_size, cur_record.name, get_output_path(outname+".scan-ED.wig"))

        elapsed_time = round((time.time() - start_time), 2)
        log_message(f"ScanFold completed. Elapsed time: {elapsed_time}s")

        # Create results README
        results_readme = os.path.join(output_folder, "RESULTS_README.md")
        with open(results_readme, 'w') as readme:
            readme.write(f"# ScanFold2.0 Results\n\n")
            readme.write(f"**Analysis completed:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            readme.write(f"## Input\n")
            readme.write(f"- **File:** {os.path.basename(myfasta)}\n")
            readme.write(f"- **Sequence:** {read_name}\n")
            readme.write(f"- **Length:** {cur_record_length} nt\n\n")
            readme.write(f"## Parameters\n")
            readme.write(f"- **Temperature:** {temperature}°C\n")
            readme.write(f"- **Window size:** {window_size} nt\n")
            readme.write(f"- **Step size:** {step_size} nt\n")
            readme.write(f"- **Shuffle type:** {shuffle_type} ({'dinucleotide' if shuffle_type == 'di' else 'mononucleotide'})\n")
            readme.write(f"- **Minimum Z-score:** {minz}\n\n")
            readme.write(f"## Output Files\n\n")
            readme.write(f"### Main Results\n")
            readme.write(f"- `{outname}` - Tab-delimited scan results with Z-scores\n")
            readme.write(f"- `{outname}.ScanFold.FinalPartners.txt` - Final base pair partners\n")
            readme.write(f"- `{outname}.ScanFold.log` - Detailed analysis log\n\n")
            readme.write(f"### Structure Files\n")
            readme.write(f"- `{basename}.no_filter.ct` - CT file with all base pairs\n")
            readme.write(f"- `{basename}.minus_1.ct` - CT file filtered at Z ≤ -1\n")
            readme.write(f"- `{basename}.minus_2.ct` - CT file filtered at Z ≤ -2\n")
            readme.write(f"- `{basename}.no_filter.dbn` - Dot-bracket notation (all pairs)\n")
            readme.write(f"- `{basename}.minus_1.dbn` - Dot-bracket notation (Z ≤ -1)\n")
            readme.write(f"- `{basename}.minus_2.dbn` - Dot-bracket notation (Z ≤ -2)\n\n")
            readme.write(f"### Track Files (WIG format)\n")
            readme.write(f"- `{outname}.scan-zscores.wig` - Z-scores track\n")
            readme.write(f"- `{outname}.scan-MFE.wig` - MFE values track\n")
            readme.write(f"- `{outname}.scan-ED.wig` - Ensemble diversity track\n")
            readme.write(f"- `{outname}.Zavg.wig` - Average Z-scores for base pairs\n\n")
            readme.write(f"### Base Pair Files\n")
            readme.write(f"- `{basename}{cur_record.name}.ALL.bp` - All base pairs with coordinates and scores\n")
            readme.write(f"- `{outname}.bp` - Filtered base pairs\n")
            readme.write(f"These files contain: i-coordinate, j-coordinate, i-nucleotide, j-nucleotide, Z-score, MFE\n")
            readme.write(f"## How to View Results\n\n")
            readme.write(f"### Quick Visualization with IGV\n\n")
            readme.write(f"**EASIEST METHOD - Use the IGV Session File:**\n\n")
            readme.write(f"```bash\n")
            readme.write(f"1. Open IGV\n")
            readme.write(f"2. File → Open Session → Select 'igv_session.json' from this folder\n")
            readme.write(f"3. Everything loads automatically!\n")
            readme.write(f"```\n\n")
            readme.write(f"#### Manual Loading (if needed):\n")
            readme.write(f"1. Download IGV from https://igv.org\n")
            readme.write(f"2. File → Load Genome from File → Select your original FASTA\n")
            readme.write(f"3. File → Load from File → Select the .wig files\n\n")
            readme.write(f"## Structure Files\n\n")
            readme.write(f"To visualize individual RNA structures:\n")
            readme.write(f"- **Dot-bracket files** (.dbn): Can be pasted into web tools like Forna or RNAfold WebServer\n")
            readme.write(f"- **CT files**: Use with VARNA or other structure visualization software\n")
            readme.write(f"- **Filtered versions**: minus_1 = Z<-1, minus_2 = Z<-2 (more stringent)\n\n")

        log_message(f"Results README created: RESULTS_README.md")

        # Create IGV session file
        igv_session_file = os.path.join(output_folder, "igv_session.json")

        # Create track list for IGV
        tracks = []

        # Add z-score track (primary, in blue)
        tracks.append({
            "name": "ScanFold Z-scores",
            "path": os.path.abspath(get_output_path(outname+".scan-zscores.wig")),
            "format": "wig",
            "color": "0,0,255",  # Blue for negative values
            "altColor": "255,0,0",  # Red for positive values
            "autoScale": False,
            "min": -3,
            "max": 1,
            "height": 100
        })

        # Add MFE track
        tracks.append({
            "name": "MFE",
            "path": os.path.abspath(get_output_path(outname+".scan-MFE.wig")),
            "format": "wig",
            "color": "0,150,0",  # Green
            "autoScale": True,
            "height": 50
        })

        # Add ED track
        tracks.append({
            "name": "Ensemble Diversity",
            "path": os.path.abspath(get_output_path(outname+".scan-ED.wig")),
            "format": "wig",
            "color": "150,0,150",  # Purple
            "autoScale": True,
            "height": 50
        })

        # Add Zavg track
        tracks.append({
            "name": "Base Pair Avg Z-scores",
            "path": os.path.abspath(get_output_path(outname+".Zavg.wig")),
            "format": "wig",
            "color": "0,100,200",  # Light blue
            "autoScale": False,
            "min": -2,
            "max": 1,
            "height": 50
        })

        # Create IGV session JSON
        igv_session = {
            "version": "1",
            "reference": {
                "fastaPath": myfasta  # Use the original input FASTA as reference
            },
            "locus": f"{read_name}:1-{cur_record_length}",  # Set view to full sequence
            "tracks": tracks
        }

        # Write IGV session file
        with open(igv_session_file, 'w') as f:
            json.dump(igv_session, f, indent=2)

        log_message(f"IGV session file created: igv_session.json")
        log_message(f"To load in IGV: File → Open Session → {igv_session_file}")
        log_message(f"All results saved to: {output_folder}")
