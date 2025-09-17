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
parser.add_argument('-r', type=int, default=100,
        help='Number of randomizations for background shuffling; default = 100')
#parser.add_argument('-a', type=str, default='rnafold',
#        help='Folding algorithm to use; default = rnafold')

### Parse arguments and convert to variables
args = parser.parse_args()
myfasta = args.filename
temperature = int(args.t)
step_size = int(args.s)
window_size = int(args.w)
randomizations = int(args.r)
#algo = str(args.a)
basename = myfasta.split('.')[0]

### 5 feature Mean MFE model
# mfe_model = tf.keras.models.load_model('/Users/ryanandrews/Desktop/scripts/5variable_meanMFEmodel')
# stddev_model = tf.keras.models.load_model('/Users/ryanandrews/Desktop/scripts/5variable_stddevMFEmodel')

### 4 feature models
mfe_model = tf.keras.models.load_model('/work/LAS/wmoss-lab/scripts/ScanFold2.0-inforna/MeanMFE')
stddev_model = tf.keras.models.load_model('/work/LAS/wmoss-lab/scripts/ScanFold2.0-inforna/StdDev')

### Start main loop
with open(myfasta, 'r') as forward_fasta:
    for cur_record in SeqIO.parse(forward_fasta, "fasta"):
        df = pd.DataFrame(columns = ["Start", "End", "Temperature", "NativeMFE",
            "Z-score", "p-value", "ED", "Sequeunce", "Structure", "centroid"])
        read_name = cur_record.name
        name = read_name
        ### Create output files
        output = str(read_name+"."+str(myfasta)+".ScanFold.")
        outname = str(read_name+".win_"+str(window_size)+".stp_"+str(step_size)+".csv")
        log_total = open(outname+".ScanFold.log", 'w')
        log_win = open(outname+".ScanFold.FinalPartners.txt", 'w')

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
         ed_list) = get_frag_feature_list(seq, step_size, window_size)

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


        meanMFE_result = mfe_model.predict(four_feature_predict)
        stddev_result = stddev_model.predict(four_feature_predict)

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
        df.to_csv(outname, sep="\t")

        minz = df["Z-score"].min()
        print(minz)
        elapsed_time = round((time.time() - start_time), 2)
        print("ScanFold-Scan complete. Elapsed time: "+str(elapsed_time)+"s")

        """
        This is the beginning of "ScanFold-Fold".
        """
        bp_dict = {}

        for index, row in df.iterrows():
            sequence = list(row["Sequence"])
            structure = list(row["Structure"])
            length = window_size
            fold_i = 0
            while fold_i < length:
                #Unpaired nucleotide
                if structure[fold_i] == '.':
                    nucleotide = sequence[fold_i]
                    coordinate = (fold_i + int(row["Start"]))
                    x = NucPair(nucleotide, coordinate, nucleotide, coordinate,
                                row["Z-score"], row["NativeMFE"], row["ED"])
                    try:
                        y = bp_dict[coordinate]
                        y.append(x)
                    except:
                        bp_dict[coordinate] = []
                        y = bp_dict[coordinate]
                        y.append(x)
                    fold_i += 1

                #Paired nucleotide
                else:
                    fold_i += 1

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
                    print("Error1")

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

                lbp_coord = int(int(lbp)+int(row["Start"])-1)
                rbp_coord = int(int(rbp)+int(row["Start"])-1)
                x = NucPair(lb, lbp_coord, rb, rbp_coord, row["Z-score"], row["NativeMFE"], row["ED"])
                z = NucPair(rb, rbp_coord, lb, lbp_coord, row["Z-score"], row["NativeMFE"], row["ED"])

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
                    t = bp_dict[rbp_coord]
                    t.append(z)
                #If i-nuc not defined, define it
                except:
                    bp_dict[rbp_coord] = []
                    t = bp_dict[rbp_coord]
                    t.append(z)
                l += 2

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
        print("ScanFold-Fold step 1 of 2 completed. Elapsed time: "+str(elapsed_time)+"s")
        print("Detecting competing pairs...")
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
        print("ScanFold-Fold step 2 of 2 completed. Elapsed time: "+str(elapsed_time)+"s")
        print("Writing files")
        strand = "forward"
        write_bp(best_bps, basename+cur_record.name+".ALL.bp", start_coordinate, name, minz)
        write_ct(final_partners, str(basename)+".no_filter.ct", float(10), strand, name, 1)
        write_ct(final_partners, str(basename)+".minus_1.ct", float(-1), strand, name, 1)
        write_ct(final_partners, str(basename)+".minus_2.ct", float(-2), strand, name, 1)
        makedbn(str(basename)+".no_filter", "NoFilter")
        makedbn(str(basename)+".minus_1", "Zavg_-1")
        makedbn(str(basename)+".minus_2", "Zavg_-2")
        write_bp(final_partners, outname+".bp", 1, name, minz)
        write_wig_dict(final_partners, outname+".Zavg.wig", name, step_size)
        write_wig(mfe_list, step_size, cur_record.name, outname+".scan-MFE.wig")
        write_wig(zscore_list, step_size, cur_record.name, outname+".scan-zscores.wig")
        write_wig((list([0] * len(mfe_list))), step_size, cur_record.name, outname+".scan-pvalue.wig")
        write_wig(ed_list, step_size, cur_record.name, outname+".scan-ED.wig")

        elapsed_time = round((time.time() - start_time), 2)
        print("ScanFold completed. Elapsed time: "+str(elapsed_time)+"s")
