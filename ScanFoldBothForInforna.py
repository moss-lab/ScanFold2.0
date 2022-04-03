### Script to run ScanFold-scan and -fold for webserver
import argparse
import logging
import os
import sys
import tempfile
import time
from datetime import datetime
from datetime import timedelta
import uuid
from ScanFoldScan import main as scan_main
from ScanFoldFold import main as fold_main
from ScanFoldFunctions import random_with_N_digits
from shutil import copyfile
import glob

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    # scan arguments
    parser.add_argument('filename',  type=str,
                        help='input filename')
    parser.add_argument('-s', type=int, default=1,
                        help='Step size; default = 1')
    parser.add_argument('-w', type=int, default=120,
                        help='Window size; default = 120')
    parser.add_argument('-r', type=int, default=100,
            help='Number of randomizations for background shuffling; default = 100')
    parser.add_argument('--algorithm', type=str, default="rnafold",
            help='Folding algorithm used; rnafold, rnastructure, mxfold')

    # fold arguments
    parser.add_argument('-f', type=int, default=-2,
                        help='filter value')
    parser.add_argument('-c', type=int, default=1,
                        help='Competition (1 for disallow competition, 0 for allow; 1 by default)')
    parser.add_argument('--id', type=str, default = "UserInput",
                        help='Name or ID of sequence being analyzed. Default "UserInput"')
    parser.add_argument('--global_refold', action='store_true', default=False,
                        help='Global refold option. Refold full sequence using Zavg <-1 and <-2 base pairs')
    parser.add_argument('--folder_name',  type=str,
                        help='Name of output folder (defaults to header name or date/time)')
    parser.add_argument('--extract', type=int, default='2',
                        help='Extract structures from minus 1 or minus 2 dbn file (2 or 1); Default = 2')
    # shared arguments
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')

    # needed for webserver
    parser.add_argument('--logfile', default=sys.stdout, type=argparse.FileType('w', encoding='UTF-8'),
            help='Path to write log file to.')
    parser.add_argument('--loglevel', default="INFO", type=str,
            help='Log level.')
    parser.add_argument('--webserver', type=str,
            help='If provided, the output folder is compressed into a tar.gz file and written to the path specified by this parameter')

    ### Parse arguments and convert to variables
    args = parser.parse_args()

    # set up logging
    loglevel = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(loglevel, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logging.basicConfig(stream=args.logfile, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=loglevel)

    try:
        # create and put results in folder based on time
        if args.folder_name == None:
            cwd = os.getcwd()
            now = datetime.now() # current date and time
            date_time = now.strftime("%m-%d-%Y_%H.%M")
            date_time = date_time + "-" + str(random_with_N_digits(3))
            folder_name = str("ScanFold_run_"+date_time)
            logging.info("\nMaking output folder named:"+folder_name)
            os.mkdir(cwd+"/"+folder_name)
            copyfile(args.filename, cwd+"/"+folder_name+"/"+args.filename)
            os.chdir(cwd+"/"+folder_name)


        if args.folder_name != None:
            folder_name = args.folder_name
            cwd = os.getcwd()
            now = datetime.now() # current date and time
            logging.info("\nMaking output folder named:"+folder_name)
            os.mkdir(cwd+"/"+folder_name)
            copyfile(args.filename, folder_name+"/"+args.filename)
            os.chdir(cwd+"/"+folder_name)
            folder_name = str(folder_name)

        # # make a temporary file to store the scan results in
        # scan_out_file = tempfile.NamedTemporaryFile(delete=False)
        # scan_out_file_name = scan_out_file.name
        # scan_out_file.close()
        #
        # fold_out_file_name = args.webserver

        # args.webserver = scan_out_file_name

        scan_main(args)

        # args.webserver = fold_out_file_name
        scan_out_file_name = glob.glob('./*.tsv')
        print(scan_out_file_name)
        args.filename = str(glob.glob('./*.tsv')[0])

        fold_main(args)

    except Exception as e:

        if args.webserver:
            # log so it shows up
            logging.error(e, exc_info=True)

        # still raise exception
        raise
