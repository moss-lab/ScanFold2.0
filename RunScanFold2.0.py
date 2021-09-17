### Script to run ScanFold-scan and -fold for dockerfile
import os
import argparse

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
parser.add_argument('--algorithm', type=str, default="rnafold",
        help='Folding algorithm used; rnafold, rnastructure, mxfold')

### ScanFold-Fold
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

args = parser.parse_args()

myfasta = args.filename
temperature = int(args.t)
step_size = int(args.s)
window_size = int(args.w)
randomizations = int(args.r)
algo = str(args.algorithm)

filter = int(args.f)
competition = int(args.c)
name = args.id
global_refold = args.global_refold
folder_name = args.folder_name
temperature = int(args.t)
algo = "rnafold"
type = "mono"


os.system('python ScanFold-Scan2.0.py {myfasta} -s {step_size} -w {window_size} -r {randomizations} -t {temperature}')
