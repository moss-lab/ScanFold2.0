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
$ python3.6 ScanFold-Scan2.0.py fasta_filename [options]

"""
import time

from ScanFoldFunctions import *
import argparse
from Bio import SeqIO
import pandas as pd
import tensorflow as tf
import os
import sys
import tf_keras as k3

import logging

def main(args):

    start_time = time.time()

    script_dir = os.path.dirname(os.path.realpath(__file__))

    myfastas = args.filename
    temperature = int(args.t)
    step_size = int(args.step)
    window_size = int(args.window)
    randomizations = int(args.r)
    algo = str(args.algorithm)

    ### 5 feature Mean MFE model
    # mfe_model = tf.keras.models.load_model('/Users/ryanandrews/Desktop/scripts/5variable_meanMFEmodel')
    # stddev_model = tf.keras.models.load_model('/Users/ryanandrews/Desktop/scripts/5variable_stddevMFEmodel')

    ### 4 feature models
    mfe_model = k3.models.load_model(os.path.join(script_dir, 'MeanMFE'))
    stddev_model = k3.models.load_model(os.path.join(script_dir, 'StdDev'))
    for myfasta in myfastas:
        basename = myfasta.split('.')[0]
        ### Start main loop
        with open(myfasta, 'r') as forward_fasta:
            for cur_record in SeqIO.parse(forward_fasta, "fasta"):
                ### Initialize dataframe
                df = pd.DataFrame(columns = ["Start", "End", "Temperature", "NativeMFE",
                    "Z-score", "p-value", "ED", "Sequence", "Structure", "centroid"])

                bp_dict = pd.DataFrame(columns = ["nuc_i", "i_bp", "nuc_j", "j_bp",
                            "z", "mfe", "ed"])

                read_name = cur_record.name
                name = read_name
                ### Create output files
                output = str(read_name+"."+str(myfasta)+".ScanFold.")
                if args.webserver:
                    outname = args.webserver
                else:
                    outname = str(read_name+".win_"+str(window_size)+".stp_"+str(step_size)+".tsv")

                ### Grab sequence and its features
                cur_record_length = len(cur_record.seq)
                seq = cur_record.seq
                seq = seq.transcribe()
                seq = seq.upper()
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
                    ed_list, di_freq_list) = get_frag_feature_list(seq, step_size, window_size, algo, temperature)

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
                df.to_csv(outname, sep="\t", index=False)

                elapsed_time = round((time.time() - start_time), 2)
                logging.info("ScanFold-Scan complete. Elapsed time: "+str(elapsed_time)+"s")

                if args.webserver:
                    break

if __name__ == "__main__":

    ### Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',  nargs="+",
                        help='input filename')
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')
    parser.add_argument('-s', '--step', type=int, default=1,
                        help='Step size; default = 1')
    parser.add_argument('-w', '--window', type=int, default=120,
                        help='Window size; default = 120')
    parser.add_argument('-r', type=int, default=100,
            help='Number of randomizations for background shuffling; default = 100')
    parser.add_argument('--algorithm', type=str, default="rnafold",
            help='Folding algorithm used; rnafold, rnastructure, mxfold')

    # needed for webserver
    parser.add_argument('--logfile', default=sys.stdout, type=argparse.FileType('w', encoding='UTF-8'),
            help='Path to write log file to.')
    parser.add_argument('--loglevel', default="INFO", type=str,
            help='Log level.')
    parser.add_argument('--webserver', type=str,
            help='If provided, only the first sequence in the fasta file is evaluated and the result is written to the path specified by this parameter')

    ### Parse arguments and convert to variables
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
