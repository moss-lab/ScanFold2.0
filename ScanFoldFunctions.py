from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import subprocess
import random
import re
from multiprocessing import get_context
import os
import statistics
import sys
sys.path.append('/work/LAS/wmoss-lab/ViennaRNA/lib/python3.6/site-packages/')
sys.path.append('/usr/local/lib/python3.6/site-packages')
sys.path.append('/usr/local/lib/python3.9/site-packages/')
import RNA
#import RNAstructure
import multiprocessing
import numpy as np
from random import randint
import tarfile
### Tensorflow imports
import pandas as pd

class NucZscore:
    #Nucleotide class; defines a nucleotide with a coordinate and a A,T,G,C,U
    def __init__(self, nucleotide, coordinate):
        self.nucleotide = nucleotide
        self.coordinate = coordinate

    def add_zscore(self, zscore):
        self.zscores.append(zscore)

    def add_pair(self, pair):
        self.pair.append(pair)

class NucPair:
    #Class to define a base pair
    def __init__(self, inucleotide, icoordinate, jnucleotide, jcoordinate, zscore, mfe, ed):
        self.inucleotide = inucleotide
        self.icoordinate = icoordinate
        self.jnucleotide = jnucleotide
        self.jcoordinate = jcoordinate
        self.zscore = zscore
        self.mfe = mfe
        self.ed = ed

class NucStructure:
    #Class to define a base pair
    def __init__(self, bond_order, coordinate, nucleotide, structure):
        self.bond_order = bond_order
        self.coordinate = coordinate
        self.nucleotide = nucleotide
        self.structure = structure

class NucStructureCount:
    #Class to define a base pair
    def __init__(self, structure_count, coordinate, nucleotide, structure):
        self.structure_count = structure_count
        self.coordinate = coordinate
        self.nucleotide = nucleotide
        self.structure = structure

class ExtractedStructure:
    def __init__(self, structure_count, sequence, structure, i, j):
        self.structure_count = structure_count
        self.sequence = sequence
        self.structure = structure
        self.i = i
        self.j = j
    def describeMe(self):
        print("object id: " + str(id(self)))
        print("structure count: " + str(self.structure_count))
        print("sequence: " + str(self.sequence))
        print("structure: " + str(self.structure))
        print("i: " + str(self.i))
        print("j: " + str(self.j))


def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def makedbn(ctfile, name):
    icoord = ()
    jcoord = ()
    kcoord = ()
    lcoord = ()
    ctfullname = ctfile+".ct"
    dbnfullname = ctfile+".dbn"
    sequence = ""
    with open(ctfullname,'r') as ct1:
        dot = ""
        data = ct1.readlines()[1:]
        for line in data:
            rows = line.split()
            icoord = int(rows[0])
            jcoord = int(rows[-2])
            #print(line)
            if len(rows) > 6:
                print("skipping header")
                continue

            elif int(rows[-2]) == 0:
                dot += '.'
                sequence += rows[1]
            else:
                sequence += rows[1]
                if icoord < jcoord:
                    with open(ctfullname,'r') as ct2:
                        data = ct2.readlines()
                        for line in data[icoord:]:
                            rows = line.split()
                            kcoord = int(rows[0])
                            lcoord = int(rows[-2])

                            if kcoord == jcoord:
                                #print(kcoord, icoord, "(")
                                dot += '('
                                break
                            else:
                                if lcoord == 0:
                                    pass
                                elif lcoord < icoord:
                                    #print('Non-nested pair found: '+str(lcoord)+' to '+str(kcoord))
                                    dot += '<'
                                    break
                                else:
                                    pass
                elif icoord > jcoord:
                    with open(ctfullname,'r') as ct2:
                        data = ct2.readlines()
                        for line in data[jcoord:]:
                            rows = line.split()
                            kcoord = int(rows[0])
                            lcoord = int(rows[-2])
                            if kcoord == icoord:
                                #print(kcoord, icoord, ")")
                                dot += ')'
                                break
                            else:
                                if lcoord == 0:
                                    pass
                                elif lcoord < jcoord:
                                    dot += '>'
                                    break
                                else:
                                    pass
                else:
                    print('Error in non-nested pair search'+'\n'+'i = '+str(icoord)+'\n'+'j = '+str(jcoord)+'\n'+'k = '+str(kcoord)+'\n'+'l = '+str(lcoord))
    #print(sequence)
    #print(dot)
    with open(dbnfullname,'w') as dbn:
        dbn.writelines(f">{name}\n")
        dbn.writelines(f"{sequence}\n")
        dbn.writelines(f"{dot}\n")

def multiprocessing(func, args,
                    workers):
    with ProcessPoolExecutor(workers) as ex:
        res = ex.map(func, args)
    return list(res)

#### Defining Dinucleotide function #####
# Taken from
# altschulEriksonDinuclShuffle.py
# P. Clote, Oct 2003
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only nonoverlapping words, presumably in a left to
# right fashion.
def computeCountAndLists(s):
  #WARNING: Use of function count(s,'UU') returns 1 on word UUU
  #since it apparently counts only nonoverlapping words UU
  #For this reason, we work with the indices.

  #Initialize lists and mono- and dinucleotide dictionaries
  List = {} #List is a dictionary of lists
  List['A'] = []; List['C'] = [];
  List['G'] = []; List['U'] = [];
  nuclList   = ["A","C","G","U"]
  s       = s.upper()
  s       = s.replace("T","U")
  nuclCnt    = {}  #empty dictionary
  dinuclCnt  = {}  #empty dictionary
  for x in nuclList:
    nuclCnt[x]=0
    dinuclCnt[x]={}
    for y in nuclList:
      dinuclCnt[x][y]=0

  #Compute count and lists
  nuclCnt[s[0]] = 1
  nuclTotal     = 1
  dinuclTotal   = 0
  for i in range(len(s)-1):
    x = s[i]; y = s[i+1]
    List[x].append( y )
    nuclCnt[y] += 1; nuclTotal  += 1
    dinuclCnt[x][y] += 1; dinuclTotal += 1
  assert (nuclTotal==len(s))
  assert (dinuclTotal==len(s)-1)
  return nuclCnt,dinuclCnt,List

def chooseEdge(x,dinuclCnt):
  numInList = 0
  for y in ['A','C','G','U']:
    numInList += dinuclCnt[x][y]
  z = random.random()
  denom=dinuclCnt[x]['A']+dinuclCnt[x]['C']+dinuclCnt[x]['G']+dinuclCnt[x]['U']
  numerator = dinuclCnt[x]['A']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['A'] -= 1
    return 'A'
  numerator += dinuclCnt[x]['C']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['C'] -= 1
    return 'C'
  numerator += dinuclCnt[x]['G']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['G'] -= 1
    return 'G'
  dinuclCnt[x]['U'] -= 1
  return 'U'

def connectedToLast(edgeList,nuclList,lastCh):
  D = {}
  for x in nuclList: D[x]=0
  for edge in edgeList:
    a = edge[0]; b = edge[1]
    if b==lastCh: D[a]=1
  for i in range(2):
    for edge in edgeList:
      a = edge[0]; b = edge[1]
      if D[b]==1: D[a]=1
  ok = 0
  for x in nuclList:
    if x!=lastCh and D[x]==0: return 0
  return 1

def eulerian(s):
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)
  #compute nucleotides appearing in s
  nuclList = []
  for x in ["A","C","G","U"]:
    if x in s: nuclList.append(x)
  #compute numInList[x] = number of dinucleotides beginning with x
  numInList = {}
  for x in nuclList:
    numInList[x]=0
    for y in nuclList:
      numInList[x] += dinuclCnt[x][y]
  #create dinucleotide shuffle L
  firstCh = s[0]  #start with first letter of s
  lastCh  = s[-1]
  edgeList = []
  for x in nuclList:
    if x!= lastCh: edgeList.append( [x,chooseEdge(x,dinuclCnt)] )
  ok = connectedToLast(edgeList,nuclList,lastCh)
  return ok,edgeList,nuclList,lastCh

def shuffleEdgeList(L):
  n = len(L); barrier = n
  for i in range(n-1):
    z = int(random.random() * barrier)
    tmp = L[z]
    L[z]= L[barrier-1]
    L[barrier-1] = tmp
    barrier -= 1
  return L

def dinuclShuffle(s):
  ok = 0
  while not ok:
    ok,edgeList,nuclList,lastCh = eulerian(s)
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)

  #remove last edges from each vertex list, shuffle, then add back
  #the removed edges at end of vertex lists.
  for [x,y] in edgeList: List[x].remove(y)
  for x in nuclList: shuffleEdgeList(List[x])
  for [x,y] in edgeList: List[x].append(y)

  #construct the eulerian path
  L = [s[0]]; prevCh = s[0]
  for i in range(len(s)-2):
    ch = List[prevCh][0]
    L.append( ch )
    del List[prevCh][0]
    prevCh = ch
  L.append(s[-1])
 # print(L)
  t = "".join(L)
  return t

#### Defining my functions #####
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'A': 'U'}
    return ''.join([complement[base] for base in dna[::-1]])

def NucleotideDictionary (lines):
    """
    Function to generate nucleotide dictionary where each key is the i
    coordinate of the nucleotide of the input sequence, and each value is a
    NucZscore class object (which contains the coordinate and nucleotide
    informations)
    """
    nuc_dict = {}
    for row in lines:
        if not row.strip():
            continue
        else:
            i = 1
            try:
                data = row.split('\t')
                icoordinate = data[0]
                sequence = simple_transcribe(str(data[7]))

            except:
                data = row.split(',')
                strand = int(data[11])
                icoordinate = data[0]
                if "A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u" in str(data[8]):
                    sequence_raw = transcribe(str(data[8]))
                elif "A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u" in str(data[7]):
                    sequence_raw = transcribe(str(data[7]))
                else:
                    raise("Could not find sequence for window")

                if strand == -1:
                    sequence = sequence_raw[::-1]

                else:
                    sequence = sequence_raw

            for nuc in sequence:
                x = NucZscore(nuc,(int(icoordinate)+int(i)-1))
                nuc_dict[x.coordinate] = x
                i += 1

    return nuc_dict;

def competing_pairs(bp_dict, coordinate):
    #Function to determine other i-nuc which compete for the same j-nuc
    comp_pairs = {}
    i = 0
    for k, v in bp_dict.items():
        if ((int(v.jcoordinate) == int(coordinate)) or
            (int(v.icoordinate) == int(coordinate))):
            x = NucPair(v.inucleotide, v.icoordinate, v.jnucleotide,
                        v.jcoordinate, v.zscore, v.mfe, v.ed)
            comp_pairs[i] = x
            i += 1
        else:
            continue

    return comp_pairs;

def best_basepair(bp_dict, nucleotide, coordinate, type):
    #Function to define best i-j pair for i-nucleotide
    zscore_dict = {}
    pair_dict = {}
    partner_key = 0
    for k, pair in sorted(bp_dict.items()):
        if int(pair.icoordinate) < int(pair.jcoordinate):
            #print("148")
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)
            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                z = pair_dict[partner_key]
                z.append(x)

            except:
                zscore_dict[partner_key] = []
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                pair_dict[partner_key] = []
                z = pair_dict[partner_key]
                z.append(x)

            sum_z = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = sum(v1)
                test = sum_z[k1] = sum(v1)

            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = statistics.mean(v1)
                test = mean_z[k1] = statistics.mean(v1)

            partner_key += 1

        elif int(pair.icoordinate) > int(pair.jcoordinate):
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)

            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                z = pair_dict[partner_key]
                z.append(x)

            except:
                zscore_dict[partner_key] = []
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                pair_dict[partner_key] = []
                z = pair_dict[partner_key]
                z.append(x)

            sum_z = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = sum(v1)
                test = sum_z[k1] = sum(v1)

            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = statistics.mean(v1)
                test = mean_z[k1] = statistics.mean(v1)

            partner_key += 1

        elif int(pair.icoordinate) == int(pair.jcoordinate):
            #print("210")
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)
            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                z = pair_dict[partner_key]
                z.append(x)

            except:
                zscore_dict[partner_key] = []
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                pair_dict[partner_key] = []
                z = pair_dict[partner_key]
                z.append(x)

            sum_z = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = sum(v1)
                test = sum_z[k1] = sum(v1)

            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = statistics.mean(v1)
                test = mean_z[k1] = statistics.mean(v1)

            partner_key += 1

        else:
            print("FAIL")
            best_bp = NucPair(pair.inucleotide, pair.icoordinate,
                              pair.jnucleotide, pair.jcoordinate, pair.zscore)

            partner_key += 1

        if type == 'sum':
            best_bp_key = min(sum_z, key = sum_z.get)
        if type == 'mean':
            best_bp_key = min(mean_z, key = mean_z.get)

    try:
        v = pair_dict[best_bp_key]
        best_bp = v[0]
    except:
        print("ERROR")
        print(k)

    return best_bp;

def write_ct(base_pair_dictionary, filename, filter, strand, name, start_coordinate):
    #Function to write connectivity table files from a list of best i-j pairs
    w = open(filename, 'w')
    w.write((str(len(base_pair_dictionary))+"\t"+name+"\n"))
    if strand == 1:
        for k, v in base_pair_dictionary.items():
            #print(start_coordinate)
            #print(v.icoordinate)
            icoordinate = str(int(v.icoordinate)-int(int(start_coordinate)-1))
            #print(icoordinate)
            jcoordinate = str(int(v.jcoordinate)-int(int(start_coordinate)-1))
            #print(jcoordinate)
            key_coordinate = str(int(k)-int(start_coordinate)+1)
            #print(key_coordinate)
            if float(v.zscore) < filter:
                if ((int(icoordinate) < int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif ((int(icoordinate) > int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif (int(icoordinate) < int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif (int(icoordinate) > int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif int(icoordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                #
                # elif (int(key_coordinate) != icoordinate) and (int(key_coordinate) != int(jcoordinate)):
                #     continue
                #     #w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("Error at", int(key_coordinate), v.inucleotide, icoordinate, v.jnucleotide, int(jcoordinate), v.zscore)
            else:
                if int(key_coordinate) == int(icoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                elif int(key_coordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    raise ValueError("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key_coordinateey_coordinateey)")
                continue

    if strand == -1:
        for k, v in sorted(base_pair_dictionary.items(), key=lambda x:x[0], reverse = True):
            # print(start_coordinate)
            # print(end_coordinate)
            # print("i="+str(v.icoordinate))
            # print("j="+str(v.jcoordinate))
            # print("k="+str(k))
            icoordinate = str(int(end_coordinate)+1-(int(int(v.icoordinate))))
            # print("i_after"+str(icoordinate))
            jcoordinate = str(int(end_coordinate)+1-(int(int(v.jcoordinate))))
            # print("j_after="+str(jcoordinate))
            key_coordinate = str(int(end_coordinate)-int(k)+1)
            # print("key="+str(key_coordinate))
            if float(v.zscore) < filter:
                if ((int(icoordinate) < int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif ((int(icoordinate) > int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif (int(icoordinate) < int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif (int(icoordinate) > int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif int(icoordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                #
                # elif (int(key_coordinate) != icoordinate) and (int(key_coordinate) != int(jcoordinate)):
                #     continue
                #     #w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("Error at", int(key_coordinate), v.inucleotide, icoordinate, v.jnucleotide, int(jcoordinate), v.zscore)
            else:
                if int(key_coordinate) == int(icoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                elif int(key_coordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    raise ValueError("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key_coordinateey_coordinateey)")
                continue

def write_dp(base_pair_dictionary, filename, filter, minz):
    #this function will create a dp file for IGV
    w = open(filename, 'w')
    for k, v in base_pair_dictionary.items():
        if float(v.zscore) < filter:
            #probability = (v.zscore/minz)
            if int(v.icoordinate) < int(v.jcoordinate):
                #w.write("%d\t%d\t%f\n" % (k, int(v.jcoordinate), float(-(math.log10(probability)))))
                w.write("%d\t%d\t%f\n" % (v.icoordinate, int(v.jcoordinate), float(((-1/minz)*(v.zscore)))/minz))
            elif int(v.icoordinate) > int(v.jcoordinate):
                w.write("%d\t%d\t%f\n" % (int(v.icoordinate), int(v.jcoordinate), float(((-1/minz)*(v.zscore)))/minz))
            elif int(v.icoordinate) == int(v.jcoordinate):
                w.write("%d\t%d\t%f\n" % (k, int(v.jcoordinate), float(((-1/minz)*(v.zscore)))/minz))
            else:
                print("Error at:", k)

def simple_transcribe(seq):
    #Function to covert T nucleotides to U nucleotides
    for ch in seq:
        rna_seq = seq.replace('T', 'U')
        return(rna_seq)### Use transcribe function via biopython (on .seq type)

def reverse_transcribe(seq):
    #Function to covert T nucleotides to U nucleotides
    for ch in seq:
        rna_seq = seq.replace('U', 'T')
        return(rna_seq)

def flip_structure(structure):
    #Function to reverse structure in a given window, for negative strand genes
    flip = {'(':')', ')':'(', '.':'.'}
    return ''.join([flip[pair] for pair in structure[::-1]])

def write_fasta(nucleotide_dictionary, outputfilename, name):
    w = open(outputfilename, 'w')
    fasta_sequence = str()
    for k, v in nucleotide_dictionary.items():
        nucleotide = v.nucleotide
        fasta_sequence += nucleotide

    w.write(">"+name+"\n")
    w.write(str(fasta_sequence)+"\n")

def write_fai (nucleotide_dictionary, filename, name):
    w = open(filename, 'w')
    name = str(name)
    length = str(len(nucleotide_dictionary))
    offset = str(utf8len(str(">"+name+"\n")))
    linebases = str(len(nucleotide_dictionary))
    linewidth = str(len(nucleotide_dictionary)+1)
    w.write("%s\t%s\t%s\t%s\t%s\n" % (name, length, offset, linebases, linewidth))

def nuc_dict_to_seq(nucleotide_dictionary):
    fasta_sequence = str()
    for k, v in nucleotide_dictionary.items():
        nucleotide = v.nucleotide
        #print(nucleotide)
        fasta_sequence += nucleotide

    return fasta_sequence

def write_wig_dict(nucleotide_dictionary, outputfilename, name, step_size, metric):
    w = open(outputfilename, 'w')
    #write wig file header
    w.write("%s %s %s %s %s\n" % ("fixedStep", "chrom="+name, "start=1", "step="+str(step_size), "span="+str(step_size)))
    #write values of zscores
    for k, v in nucleotide_dictionary.items():
        if str(metric) == 'zscore':
            zscore = float(v.zscore)
            w.write("%f\n" % (zscore))

        elif str(metric) == 'mfe':
            mfe = float(v.mfe)
            w.write("%f\n" % (mfe))

        elif str(metric) == 'ed':
            ed = float(v.ed)
            w.write("%f\n" % (ed))
        else:
            print("Set metric to zscore, mfe, or ed")

def write_wig(metric_list, step, name, outputfilename):
    w = open(outputfilename, 'w')
    w.write("%s %s %s %s %s\n" % ("fixedStep", "chrom="+name, "start=1", "step="+str(step), "span="+str(step)))
    for metric in metric_list:
        # try:
        #     float(metric)
        #     print("GOOD", metric)
        # except:
        #     print(metric)
        if metric == "#DIV/0!":
            w.write("%s\n" % (metric))
        else:
            try:
                w.write("%f\n" % (metric))
            except:
                w.write("%s\n" % (metric))
                #print(metric)

def write_bp(base_pair_dictionary, filename, start_coordinate, name, minz):

    w = open(filename, 'w')
        #set color for bp file (igv format)
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 55, 129, 255, str("Less than -2 "+str(minz))))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 89, 222, 111, "-1 to -2"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 236, 236, 136, "0 to -1"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 199, 199, 199, "0"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 228, 228, 228, "0 to 1"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 243, 243, 243, "1 to 2"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 247, 247, 247, str("Greater than 2")))

    i = 0
    for k, v in base_pair_dictionary.items():
        #choose color
        if float(v.zscore) < float(-2):
            score = str(0)
            #print(k, v.zscore, score)

        elif (float(v.zscore) < int(-1)) and (float(v.zscore) >= -2):
            score = str(1)
            #print(k, v.zscore, score)

        elif (float(v.zscore) < int(0)) and (float(v.zscore) >= -1):
            score = str(2)
            #print(k, v.zscore, score)

        elif float(v.zscore) == 0 :
            score = str(3)
            #print(k, v.zscore, score)

        elif 0 < float(v.zscore) <= 1:
            score = str(4)
            #print(k, v.zscore, score)

        elif 1 < float(v.zscore) <= 2:
            score = str(5)
            #print(k, v.zscore, score)

        elif float(v.zscore) > 2:
            score = str(6)
            #print(k, v.zscore, score)

        else:
            print(k, v.zscore, score)


        score = str(score)

        # ensure coordinates to start at 1 to match with converted fasta file
        sc = int(int(start_coordinate)-1)
        #print(length)


        if int(v.icoordinate) < int(v.jcoordinate):
            #w.write("%d\t%d\t%f\n" % (k, int(v.jcoordinate), float(-(math.log10(probability)))))
            # print(name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score)
            # print(name, str(int(v.icoordinate)-sc), str(int(v.icoordinate)-sc), str(int(v.jcoordinate)-sc), str(int(v.jcoordinate)-sc), score)
            w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, int(v.icoordinate)-sc, int(v.icoordinate)-sc, int(v.jcoordinate)-sc, int(v.jcoordinate)-sc, score))
        elif int(v.icoordinate) > int(v.jcoordinate):
            # print(name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score)
            # print(name, str(int(v.icoordinate)-sc), str(int(v.icoordinate)-sc), str(int(v.jcoordinate)-sc), str(int(v.jcoordinate)-sc), score)
            w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, int(v.icoordinate)-sc, int(v.icoordinate)-sc, int(v.jcoordinate)-sc, int(v.jcoordinate)-sc, score))
        elif int(v.icoordinate) == int(v.jcoordinate):
            # print(name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score)
            # print(name, str(int(v.icoordinate)-sc), str(int(v.icoordinate)-sc), str(int(v.jcoordinate)-sc), str(int(v.jcoordinate)-sc), score)
            w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, k-sc, k-sc, int(v.jcoordinate)-sc, int(v.jcoordinate)-sc, score))
        else:
            print("2 Error at:", k)

def utf8len(s):
    return len(s.encode('utf-8'))

def write_fai (nucleotide_dictionary, filename, name):
    w = open(filename, 'w')
    name = str(name)
    length = str(len(nucleotide_dictionary))
    offset = str(utf8len(str(">"+name+"\n")))
    linebases = str(len(nucleotide_dictionary))
    linewidth = str(len(nucleotide_dictionary)+1)
    w.write("%s\t%s\t%s\t%s\t%s\n" % (name, length, offset, linebases, linewidth))

###### Function to calculate ZScore on list of MFEs #################
def pvalue_function(energy_list, randomizations):
    below_native = 0
    total_count = len(energy_list)
    native_mfe = float(energy_list[0])
    #scrambled_mean_mfe = statistics.mean(energy_list[1:randomizations])
    for MFE in energy_list:
        if float(MFE) < float(native_mfe):
            below_native += 1

    pvalue = float(float(below_native) / float(total_count))

    return pvalue;

###### Function to calculate ZScore on list of MFEs #################
def zscore_function(energy_list, randomizations):
    mean = np.mean(energy_list)
    sd = np.std(energy_list)
    native_mfe = energy_list[0]
    scrambled_mean_mfe = np.mean(energy_list[1:randomizations])
    #scrambled_sd = statistics.stdev(energy_list[1:randomizations])
    if sd != 0:
        zscore = (native_mfe - scrambled_mean_mfe)/sd
    if sd == 0:
        zscore = float(00.00)
    return zscore;

def flip_structure(structure):
    #Function to reverse structure in a given window, for negative strand genes
    flip = {'(':')',  ')':'(',  '.':'.',  '&':'&'}
    return ''.join([flip[pair] for pair in structure[::-1]])

def rna_fold(frag, temperature):
    try:
        if algo == "rnastructure":
            p = RNAstructure.RNA.fromString(str(frag))
            p.FoldSingleStrand(mfeonly=True)
            MFE = p.GetFreeEnergy(1)
            #(structure, MFE) = RNA.fold(str(frag))

        if algo == "rnafold":
            fc = RNA.fold_compound(str(frag), md)
            (structure, MFE) = fc.mfe()
            #print(md.temperature)
        return MFE;

    except:

        args = ["RNAfold", "-p", "-T", str(temperature)]
        fc = subprocess.run(args, input=str(frag), check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = str(fc.stdout)
        test = out.splitlines()
        structure = test[1].split()[0]
        centroid = test[3].split()[0]
        MFE = test[1].split(" ", 1)[1]
        try:
            MFE = float(re.sub('[()]', '', MFE))
        except:
            print("Error parsing MFE values", test)
        ED = float(test[4].split()[-1])

        return (structure, centroid, MFE, ED)

def rna_refold(frag, temperature, constraint_file):
    args = ["RNAfold", "-p", "-T", str(temperature), '-C', constraint_file]
    fc = subprocess.run(args, input=frag, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = str(fc.stdout)
    test = out.splitlines()
    structure = test[2].split()[0]
    centroid = test[4].split()[0]
    MFE = test[2].split(" ", 1)[1]
    try:
        MFE = float(re.sub('[()]', '', MFE))
    except:
        print("Error parsing MFE values", test)
    ED = float(test[5].split()[-1])

    return (structure, centroid, MFE, ED)

def rna_folder(arg):
    (frag, temperature, algo) = arg
    md = RNA.md()
    md.temperature = int(temperature)
    if algo == "rnastructure":
        p = RNAstructure.RNA.fromString(str(frag))
        p.FoldSingleStrand(mfeonly=True)
        MFE = p.GetFreeEnergy(1)
        #(structure, MFE) = RNA.fold(str(frag))

    if algo == "rnafold":
        fc = RNA.fold_compound(str(frag), md)
        (structure, MFE) = fc.mfe()

    return MFE;



# def rna_cofolder(frag1, frag2):
#     frag1 = str(frag1)
#     frag2 = str(frag2)
#     (structure, energy) = RNA.duplex(frag1, frag2)
#     #print(md.temperature)
#     return energy;

def randomizer(frag):
    result = ''.join(random.sample(frag,len(frag)))
    return result;

###### Function to calculate MFEs using RNAfold #################
def energies(seq_list, temperature, algo):
    energy_list = []

    try:
        energy_list = multiprocessing(rna_folder, [(sequence, temperature, algo) for sequence in seq_list], 12)
    except:
        for sequence in seq_list:
            #fc = RNA.fold_compound(str(sequence))
            (structure, MFE) = RNA.fold(str(sequence)) # calculate and define variables for mfe and structure
            energy_list.append(MFE) # adds the native fragment to list

    return energy_list;

###### Function to calculate MFEs using RNAcofold #################
def cofold_energies(frag1, seq_list):
    energy_list = []
    for sequence in seq_list:
        scrambled_frag1 = ''.join(random.sample(frag1,len(frag1)))
        #scrambled_frag1 = frag1
        #print(scrambled_frag1)
        #fc = RNA.fold_compound(str(sequence))
        duplex = RNA.duplexfold(str(scrambled_frag1), str(sequence)) # calculate and define variables for mfe and structure
        # duplex = RNA.duplexfold(str(frag1), str(sequence)) # calculate and define variables for mfe and structure
        MFE = duplex.energy
        energy_list.append(MFE) # adds the native fragment to list

    return energy_list;


######Function to create X number of scrambled RNAs in list #################
#test
def scramble(text, randomizations, type):
    frag = str(text)
    frag_seqs = []
    if type == "di":
        frag = simple_transcribe(frag)
        for _ in range(randomizations):
            result = dinuclShuffle(frag)
            frag_seqs.append(result)
    elif type == "mono":
        try:
            frag_seqs = multiprocessing(randomizer, [frag for i in range(randomizations)], 12)
        except:
            for _ in range(int(randomizations)):
                result = ''.join(random.sample(frag,len(frag)))
                frag_seqs.append(result)
    else:
        print("Shuffle type not properly designated; please input \"di\" or \"mono\"")

    return frag_seqs;

def get_structure(rnastructure_object):
    structure = []
    p = rnastructure_object
    for i in range(1, p.GetSequenceLength()+1):
        pair = p.GetPair(i)
        if pair == 0:
            structure.append(".")
        elif pair > i:
            structure.append("(")
        elif pair < i:
            structure.append(")")
        else:
            raise ValueError('"get_structure" funciton error: RNA structure class could not be read to generate folding structure')

    return structure;

def findpair(nucdict, k):
    #print(k)
    base = nucdict[k].coordinate
    order = nucdict[k].bond_order
    #print(base, nucdict[k].nucleotide)
    #print("Finding pair for ", nucdict[k].coordinate, nucdict[k].structure)
    if nucdict[k].structure == ".":
        print("Cant find pair for unpaired nucleotide at i =", k+1)
    if nucdict[k].structure == "(":
        i = 1
        test = 1
        while test > 0:

            if nucdict[k+i].structure == ".":
                i += 1

            if nucdict[k+i].structure == "(":
                i += 1
                test +=1

            if nucdict[k+i].structure == ")":
                test -= 1
                if test == 0:
                    #print("FOUND")
                    pair = nucdict[k+i]
                    #print(nucdict[k].coordinate, nucdict[k].structure, pair.structure, pair.coordinate)
                i += 1

            #print(i)

    if nucdict[k].structure == ")":
        i = 1
        test = 1
        #print(nucdict[k].coordinate)
        while test > 0:
            #print(test)
            if nucdict[k-i].structure == ".":
                i += 1

            if nucdict[k-i].structure == ")":
                i += 1
                test += 1

            if nucdict[k-i].structure == "(":
                test -= 1
                if test == 0:
                    #print("FOUND")
                    pair = nucdict[k-i]
                    #print(test, nucdict[k].coordinate, nucdict[k].structure, pair.structure, pair.coordinate)
                i += 1

    return pair;

def dbn2ct(dbnfile):
    #print("Running dbn2ct")
    #Function to write connectivity table files from a list of best i-j pairs
    with open(dbnfile, 'r') as f:
        lines = f.readlines()
        sequence_raw = str(lines[1].strip())
        structure_raw = str(lines[2].strip())
        sequence = list(sequence_raw)
        structure = list(structure_raw)
        length = len(sequence)
        length_st = len(structure)
        bond_order = []
        bond_count = 0
        nuc_dict = {}
        #print(sequence, structure, length)
        #Inititate base pair tabulation variables
        bond_order = []
        bond_count = 0
        nuc_dict = {}
        #print(length)
        if length != length_st:
            print(sequence, structure)
            raise("ERROR structure and sequence not same length.")


        #Iterate through sequence to assign nucleotides to structure type
        m = 0
        #print(length)
        while  m <= length-1:
            #print(structure[m])
            if structure[m] == '(':
                # print(m, structure[m])
                bond_count += 1
                bond_order.append(bond_count)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif structure[m] == ')':
                # print(m, structure[m])
                bond_order.append(bond_count)
                bond_count -= 1
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif str(structure[m]) == ( '.' ):
                # print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '<' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '>' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '{' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '}' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif str(structure[m]) == ( '[' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( ']' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            else:
                print("Error", bond_count, (m+1), sequence[m], structure[m])
                m += 1
                continue
                #print("no")



    # for k, v in nuc_dict.items():
    #     print(str(k+1), str(v.bond_order), str(v.nucleotide), str(k), str(k-1), str(k+1), str(0), str(k+1))
    ct_out = dbnfile.replace('.dbn', '.ct')
    with open(ct_out,'w') as ct:
        ct.write(str(length-1)+"\t"+str("SequenceID"+"\n"))
        for k, v in nuc_dict.items():
            #print(v.coordinate, v.structure)
            if v.bond_order >= 0:
                if v.structure == ".":
                    ct.write("%d %s %d %d %d %d\n" % ((k+1), v.nucleotide, k, k+2, 0, k+1))

                if v.structure == "(":
                    pair = findpair(nuc_dict, k)
                    ct.write("%d %s %d %d %d %d\n" % ((k+1), v.nucleotide, k, k+2, pair.coordinate, k+1))

                if v.structure == ")":
                    pair = findpair(nuc_dict, k)
                    ct.write("%d %s %d %d %d %d\n" % (v.coordinate, v.nucleotide, k, k+2, pair.coordinate, k+1))

def random_with_N_digits(n):
    range_start = 10**(n-1)
    range_end = (10**n)-1
    return randint(range_start, range_end)

def get_gc_content(frag):
    ###Ensure frag is string
    frag = str(frag)
    if 'C' and 'G' in frag:
      A_count = frag.count("A")+frag.count("a")
      G_count = frag.count("G")+frag.count("g")
      C_count = frag.count("C")+frag.count("c")
      T_count = frag.count("T")+frag.count("t")+frag.count("U")+frag.count("u")
      gc_content = round(float(G_count+C_count)/float(A_count+T_count+G_count+C_count),5)
    else:
      gc_content = 0

    return gc_content

def get_cg_ratio(frag):
    frag = str(frag)
    if 'C' and 'G' in frag:
      A_count = frag.count("A")+frag.count("a")
      G_count = frag.count("G")+frag.count("g")
      C_count = frag.count("C")+frag.count("c")
      T_count = frag.count("T")+frag.count("t")+frag.count("U")+frag.count("u")
      cg_ratio = C_count/(C_count+G_count)
    else:
      cg_ratio = 0

    return cg_ratio

def get_au_ratio(frag):
    frag = str(frag)
    if 'A' and 'U' in frag:
      A_count = frag.count("A")+frag.count("a")
      G_count = frag.count("G")+frag.count("g")
      C_count = frag.count("C")+frag.count("c")
      T_count = frag.count("T")+frag.count("t")+frag.count("U")+frag.count("u")
      au_ratio = A_count/(A_count+T_count)
    else:
      au_ratio = 0

    return au_ratio

def get_svm_zscore(frag):
    frag = str(frag)
    fc = RNA.fold_compound(frag)
    (structure, mfe) = fc.mfe()
    gc_content = get_gc_content(frag)
    cg_ratio = get_cg_ratio(frag)
    au_ratio = get_au_ratio(frag)
    length = len(frag)
    #node_mono = {1:gc_content, 2:cg_ratio, 3:au_ratio, 4:length}
    node_mono = ((1, gc_content), (2, cg_ratio), (3, au_ratio), (4, length))
    y_node_mono = ((1, gc_content), (2, cg_ratio), (3, au_ratio), (4, length))
    print(type(node_mono))
    avg_m = svm_load_model('/Users/ryanandrews/Desktop/scripts/RNAz/models/mfe_avg.model')
    avg = svm_predict(node_mono, node_mono, avg_m)
    stdv_m = svm_load_model('/Users/ryanandrews/Desktop/scripts/RNAz/models/mfe_stdv.model')
    stdv = svm_predict(node_mono, node_mono, stdv_m)
    print(mfe, avg, stdv)
    svm_zscore = mfe-avg/stdv

    return svm_zscore

def get_di_freqs(frag):
    ### code taken from https://pythonforbiologists.com/dictionaries
    frag = str(frag)
    frag_list = [frag[i:i+2] for i in range(0, len(frag))]
    dinucleotides = ['AA','AU','AG','AC',
                     'UA','UU','UG','UC',
                     'GA','GU','GG','GC',
                     'CA','CU','CG','CC']
    all_counts = []
    for dinucleotide in dinucleotides:
        count = frag_list.count(dinucleotide)
        #print("count is " + str(count) + " for " + dinucleotide)
        all_counts.append(count/(len(frag)-1))

    return(all_counts)

def get_dinucleotide_counts(frag):
    ### code taken from https://pythonforbiologists.com/dictionaries
    frag = str(frag)
    frag_list = [frag[i:i+2] for i in range(0, len(frag))]
    #print(frag_list)
    dinucleotides = ['AA','AU','AG','AC',
                     'UA','UU','UG','UC',
                     'GA','GU','GG','GC',
                     'CA','CU','CG','CC']
    all_counts = []
    for dinucleotide in dinucleotides:
        count = frag_list.count(dinucleotide)
        #print("count is " + str(count) + " for " + dinucleotide)
        all_counts.append(count)

    return(all_counts)

def get_frag_feature_list(seq, step_size, window_size, algo, temperature):
    frag_list = []
    start_coordinate_list = []
    end_coordinate_list = []
    GCpercent_list = []
    CGratio_list = []
    AUratio_list = []
    structure_list = []
    centroid_list = []
    mfe_list = []
    ed_list = []
    di_freq_list = []

    i = 0
    while i == 0 or i <= (len(seq) - window_size):
        start_coordinate = i+1 # This will just define the start nucleotide coordinate value
        start_coordinate_list.append(start_coordinate)
        end_coordinate = i + window_size
        end_coordinate_list.append(end_coordinate)
        frag = str(seq[i:i+int(window_size)]) # This breaks up sequence into fragments
        #print(frag)
        frag_list.append(frag)
        GCpercent = get_gc_content(frag)
        CGratio = get_cg_ratio(frag)
        AUratio = get_au_ratio(frag)
        di_freqs = get_di_freqs(frag)
        di_freq_list.append(di_freqs)
        GCpercent_list.append(GCpercent)
        CGratio_list.append(CGratio)
        AUratio_list.append(AUratio)
        if algo == "rnafold":
            md = RNA.md()
            md.temperature = temperature
            fc = RNA.fold_compound(frag,md)
            fc.pf()
            (centroid, distance) = fc.centroid()
            ensemble_diversity = round(fc.mean_bp_distance(), 2)
            (structure, native_mfe) = fc.mfe()
            #print(structure, native_mfe)
            native_mfe = round(native_mfe, 2)
        elif algo == "rnastructure":
            temp_kelvin = temperature+273.15
            args = ["Fold", "-k", "-mfe", "-T", str(temp_kelvin), "-", "-"]
            fc = subprocess.run(args, input=str(frag), check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = str(fc.stdout)
            test = out.splitlines()
            structure = test[2].split()[0]
            #print(structure)
            centroid = "NA"
            native_mfe = test[0].split(" ", 3)[2]
            try:
                native_mfe = float(native_mfe)
            except:
                print("Error parsing MFE values", test)
            ensemble_diversity = 0.0
        elif algo == "test":
            structure = "[[.[[[[]]...]]]]..[[[..((((.......]]].))))((((([[[[{<)))))....]]]]}>..............................[[[[[[[....]].]]]]]..."
            native_mfe = -20
            centroid = "NA"
            ensemble_diversity = 0.0
        else:
            print("No folding algorihm properly passed to function.")

        centroid_list.append(centroid)
        ed_list.append(ensemble_diversity)
        structure_list.append(structure)
        mfe_list.append(native_mfe)
        i += step_size
        #print(i)

    if step_size > 1:
        start_coordinate = (len(seq)-window_size)+1 # This will just define the start nucleotide coordinate value
        start_coordinate_list.append(start_coordinate)
        end_coordinate = len(seq)
        end_coordinate_list.append(end_coordinate)
        frag = str(seq[start_coordinate:end_coordinate]) # This breaks up sequence into fragments
        #print(frag)
        frag_list.append(frag)
        GCpercent = get_gc_content(frag)
        CGratio = get_cg_ratio(frag)
        AUratio = get_au_ratio(frag)
        di_freqs = get_di_freqs(frag)
        di_freq_list.append(di_freqs)
        GCpercent_list.append(GCpercent)
        CGratio_list.append(CGratio)
        AUratio_list.append(AUratio)
        fc = RNA.fold_compound(frag)
        fc.pf()
        (centroid, distance) = fc.centroid()
        ensemble_diversity = round(fc.mean_bp_distance(), 2)
        (structure, native_mfe) = fc.mfe()
        #print(structure, native_mfe)
        native_mfe = round(native_mfe, 2)
        centroid_list.append(centroid)
        ed_list.append(ensemble_diversity)
        structure_list.append(structure)
        mfe_list.append(native_mfe)

    length_list = [window_size] * len(frag_list)
    return (start_coordinate_list, end_coordinate_list, frag_list, length_list,
            GCpercent_list, CGratio_list, AUratio_list,
            mfe_list, structure_list, centroid_list, ed_list, di_freq_list)

def write_bp_from_df(dataframe, filename, start_coordinate, name, minz):

    w = open(filename, 'w')
        #set color for bp file (igv format)
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 55, 129, 255, str("Less than -2 "+str(minz))))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 89, 222, 111, "-1 to -2"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 236, 236, 136, "0 to -1"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 199, 199, 199, "0"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 228, 228, 228, "0 to 1"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 243, 243, 243, "1 to 2"))
    w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 247, 247, 247, str("Greater than 2")))

    for index, row in dataframe.iterrows():
        #choose color
        if float(row["z_avg"]) < float(-2):
            score = str(0)
            #print(k, row["z_avg"], score)

        elif (float(row["z_avg"]) < int(-1)) and (float(row["z_avg"]) >= -2):
            score = str(1)
            #print(k, row["z_avg"], score)

        elif (float(row["z_avg"]) < int(0)) and (float(row["z_avg"]) >= -1):
            score = str(2)
            #print(k, row["z_avg"], score)

        elif float(row["z_avg"]) == 0 :
            score = str(3)
            #print(k, row["z_avg"], score)

        elif 0 < float(row["z_avg"]) <= 1:
            score = str(4)
            #print(k, row["z_avg"], score)

        elif 1 < float(row["z_avg"]) <= 2:
            score = str(5)
            #print(k, row["z_avg"], score)

        elif float(row["z_avg"]) > 2:
            score = str(6)
            #print(k, row["z_avg"], score)

        else:
            print(f"Error finding zavg value at {row}")


        score = str(score)

        # ensure coordinates to start at 1 to match with converted fasta file
        sc = int(int(start_coordinate)-1)
        #print(length)


        if int(row["i_bp"]) < int(row["j_bp"]):
            #w.write("%d\t%d\t%f\n" % (k, int(row["j_bp"]), float(-(math.log10(probability)))))
            # print(name, int(row["i_bp"]), int(row["i_bp"]), int(row["j_bp"]), int(row["j_bp"]), score)
            # print(name, str(int(row["i_bp"])-sc), str(int(row["i_bp"])-sc), str(int(row["j_bp"])-sc), str(int(row["j_bp"])-sc), score)
            w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, int(row["i_bp"])-sc, int(row["i_bp"])-sc, int(row["j_bp"])-sc, int(row["j_bp"])-sc, score))
        elif int(row["i_bp"]) > int(row["j_bp"]):
            # print(name, int(row["i_bp"]), int(row["i_bp"]), int(row["j_bp"]), int(row["j_bp"]), score)
            # print(name, str(int(row["i_bp"])-sc), str(int(row["i_bp"])-sc), str(int(row["j_bp"])-sc), str(int(row["j_bp"])-sc), score)
            w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, int(row["i_bp"])-sc, int(row["i_bp"])-sc, int(row["j_bp"])-sc, int(row["j_bp"])-sc, score))
        elif int(row["i_bp"]) == int(row["j_bp"]):
            # print(name, int(row["i_bp"]), int(row["i_bp"]), int(row["j_bp"]), int(row["j_bp"]), score)
            # print(name, str(int(row["i_bp"])-sc), str(int(row["i_bp"])-sc), str(int(row["j_bp"])-sc), str(int(row["j_bp"])-sc), score)
            w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, (row["i_bp"])-sc, (row["i_bp"])-sc, int(row["j_bp"])-sc, int(row["j_bp"])-sc, score))
        else:
            print("2 Error at:", row["i_bp"])

def write_ct_from_df(dataframe, filename, filter, strand, name, start_coordinate):
    #Function to write connectivity table files from a list of best i-j pairs
    w = open(filename, 'w')
    w.write((str(len(dataframe))+"\t"+name+"\n"))
    if strand == 1:
        i = 0
        for index, row in dataframe.iterrows():
            #print(start_coordinate)
            #print(row["i_bp"])
            icoordinate = str(int(row["i_bp"])-int(int(start_coordinate)-1))
            #print(icoordinate)
            jcoordinate = str(int(row["j_bp"])-int(int(start_coordinate)-1))
            #print(jcoordinate)
            key_coordinate = (icoordinate)
            #print(key_coordinate, icoordinate, jcoordinate)
            if float(row["z_avg"]) < filter:
                if ((int(icoordinate) < int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif ((int(icoordinate) > int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif (int(icoordinate) < int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_j"], int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif (int(icoordinate) > int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_j"], int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif int(icoordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                #
                # elif (int(key_coordinate) != icoordinate) and (int(key_coordinate) != int(jcoordinate)):
                #     continue
                #     #w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("Error at", int(key_coordinate), row["nuc_i"], icoordinate, row["nuc_j"], int(jcoordinate), row["z_avg"])
                i += 1
            else:
                #print(key_coordinate, icoordinate, jcoordinate)
                if int(key_coordinate) == int(icoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                elif int(key_coordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_j"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key_coordinateey_coordinateey)")
                i += 1

    if strand == -1:
        print("FIX REVERSE STRAND WRITE CT FROM DF")
        for index, row in sorted(dataframes.iterrows(), key=lambda x:x[0], reverse = True):
            # print(start_coordinate)
            # print(end_coordinate)
            # print("i="+str(v.icoordinate))
            # print("j="+str(row["j_bp"]))
            # print("k="+str(k))
            icoordinate = str(int(end_coordinate)+1-(int(int(row["i_bp"]))))
            # print("i_after"+str(icoordinate))
            jcoordinate = str(int(end_coordinate)+1-(int(int(row["j_bp"]))))
            # print("j_after="+str(jcoordinate))
            key_coordinate = str(int(end_coordinate)-int(k)+1)
            # print("key="+str(key_coordinate))
            if float(row["z_avg"]) < filter:
                if ((int(icoordinate) < int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif ((int(icoordinate) > int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif (int(icoordinate) < int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_j"], int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif (int(icoordinate) > int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_j"], int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif int(icoordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                #
                # elif (int(key_coordinate) != icoordinate) and (int(key_coordinate) != int(jcoordinate)):
                #     continue
                #     #w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("Error at", int(key_coordinate), row["nuc_i"], icoordinate, row["nuc_j"], int(jcoordinate), row["z_avg"])
            else:
                if int(key_coordinate) == int(icoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_i"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                elif int(key_coordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), row["nuc_j"], int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    raise ValueError("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key_coordinateey_coordinateey)")
                continue

def merge_files(destination, *sources):
    with open(destination, 'w') as outfile:
        for fname in sources:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def make_tar(destination, source):
    with open(destination, 'wb') as output_wb:
        tar = tarfile.open(mode="w:gz", fileobj=output_wb)
        _, _, filenames = next(os.walk(source))
        for name in filenames:
            tar.add(name)
        tar.close()

def random_with_N_digits(n):
    range_start = 10**(n-1)
    range_end = (10**n)-1
    return randint(range_start, range_end)

def structure_extract(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=int, default=100,
                        help='randomizations for z-score calculation')
    parser.add_argument('-f', type=int, default=0,
                        help='number of flanking nucleotides')
    parser.add_argument('--type', type=str, default='mono',
                        help='Randomization type')
    parser.add_argument('--filename',  type=str, default='Zavg_-2_pairs.dbn',
                        help='input filename')
    parser.add_argument('--structure_extract_file', type=str, default = "ExtractedStructures.gff3",
                        help='structure_extract_file path')
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')
    parser.add_argument('--algo', type=str, default='rnafold',
                        help='Select RNA folding algorithm')
    parser.add_argument('--name', type=str,
                        help='Name or ID of sequence being analyzed. Default "UserInput"')
    args = parser.parse_args()
    filename = args.filename
    randomizations = int(args.r)
    temperature = int(args.t)
    flanking = int(args.f)
    structure_extract_file = args.structure_extract_file
    algo = str(args.algo)
    type = str(args.type)
    name = args.name

    with open(filename, 'r') as f:
        lines = f.readlines()
        header = lines[0]
        split_header = header.split('\t')
        #print(split_header[10])
        accesion = filename
        #print (accesion)
        fname_list = accesion.split('/')
        fname_part = fname_list[len(fname_list) -1]
        #print(fname_part)
        fname = f'ExtrStr_{fname_part}.txt'
        MAIN = os.path.join(root, fname)

        se = open(fname, 'w')
        se.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("StructureNumber", "i", "j", "sequence", "ScanFoldPairs", "RefoldedStructure", "RefoldedMFE", "Refoldedz-score", "RefoldedED"))

        #print(title)
        sequence_raw = str(lines[1])
        #print(sequence_raw)
        structure_raw = str(lines[2])
        #print(structure_raw)


        bp_dict = {}

        sequence = list(sequence_raw)
        structure = list(structure_raw)

        length = len(sequence)
        length_st = len(structure)
        #print(length, length_st)

        #Inititate base pair tabulation variables
        bond_order = []
        bond_count = 0
        nuc_dict = {}
        #Iterate through sequence to assign nucleotides to structure type
        m = 0
        #print(length)
        while  m < length-1:
            #print(m)
            if structure[m] == '(':
                #print(m, structure[m])
                bond_count += 1
                bond_order.append(bond_count)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif structure[m] == ')':
            #    print(m, structure[m])
                bond_order.append(bond_count)
                bond_count -= 1
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif str(structure[m]) == ( '.' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '<' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '>' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '{' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '}' ):
            #    print(m, structure[m])
                bond_order.append(0)
                nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            else:
                #print("Error", bond_count, (m+1), sequence[m], structure[m])
                m += 1
                continue
                # print("no")
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

            elif str(structure[m]) == ( '[' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( ']' ):
            #    print(m, structure[m])
                bond_order_pk.append(0)
                nuc_dict_pk[m] = NucStructure(bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
            else:
                #print("Error", bond_count_pk, (m+1), sequence[m], structure[m])
                m += 1
                continue

        #print(bond_order)
        #Initiate base_pair list
        base_pairs = []

        #Create empty variable named test
        test = ""

        #Iterate through bond order
        j = 0
        structure_count = 0
        structure_end = []
        structure_start = []
        while j < length:
            #print(f"{j} nested")
            try:
                if (nuc_dict[j].bond_order == 1) and (nuc_dict[j].structure == '('):
                    structure_count += 1
                    #print(nuc_dict[j].structure, j)
                    structure_start.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict[j].nucleotide, nuc_dict[j].structure))
                    j += 1

                elif (nuc_dict[j].bond_order == 0) and (nuc_dict[j].structure == ')'):
                    structure_count += 1
                    #print(nuc_dict[j].structure, j)
                    structure_end.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict[j].nucleotide, nuc_dict[j].structure))
                    j += 1
                else:
                    j += 1
            except:
                j += 1
                continue


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

        #print(structure_start[0].coordinate, structure_end[0].coordinate)

        extracted_structure_list = []


        l = 0
        while l < int(len(structure_start)):
            offset = flanking
            s = structure_start_coordinate =  int((structure_start[l].coordinate)-offset-1)
            e = structure_end_coordinate = int((structure_end[l].coordinate)+offset-1)

            seq = ""
            fold = ""
            for k, v in nuc_dict.items():
                if s <= k <= e:
                    seq += str(v.nucleotide)
                    fold += str(v.structure)

            extracted_structure_list.append(ExtractedStructure(l, seq, fold, s, e))

            l += 1

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
        #try:
            #print(str(i))
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
            #print(structure)
            #fmfe = fc.pbacktrack()
            #print(str(fmfe))
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
            #print(pvalue)
            pvalue_total.append(pvalue)
            accession = str(name)
            ps_title = f"motif_{motif_num} coordinates {es.i} - {es.j}"
            #print(macro)

            with open(f"{name}_motif_{motif_num}.dbn", 'w') as es_dbn:
                es_dbn.write(f">{name}_motif_{motif_num}_coordinates:{es.i}-{es.j}\n{frag}\n{MFE_structure}")
            dbn2ct(f"{name}_motif_{motif_num}.dbn")
            ###Create figure using forgi
            # with open('tmp.fa', 'w') as file:
            #     file.write('>tmp\n%s\n%s' % (frag, MFE_structure))
            # cg = forgi.load_rna("tmp.fa", allow_many=False)
            # fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7, backbone_kwargs={"linewidth":3})
            # plt.savefig(ps_title+".png")
            # plt.clf()


            RNA.PS_rna_plot_a(frag, MFE_structure, "motif_"+str(motif_num)+".ps", '', '')
            gff_attributes = f'motif_{motif_num};sequence={es.sequence};structure={str(es.structure)};refoldedMFE={str(MFE_structure)};MFE(kcal/mol)={str(MFE)};z-score={str(zscore)};ED={str(ED)}'
            se.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str("."), str("RNA_sequence_secondary_structure"), str(int((es.i)+1)), str(int((es.j)+1)), str("."), str("."), str("."), gff_attributes))
            motif_num += 1

        se.close()

# def pairs_to_dataframe(sequence, structure):
#     fold_i = 0
#     for pair in structure:
#         #Unpaired nucleotide
#         if pair == '.':
#             nucleotide = sequence[fold_i]
#             coordinate = (fold_i + start_coordinate)
#             #print(nucleotide, coordinate)
#             '''
#             bp_dict = bp_dict.append({'nuc_i':nucleotide, 'i_bp':coordinate,
#                                 'nuc_j':nucleotide, 'j_bp':coordinate,
#                         'z':row["Z-score"], 'mfe':row["NativeMFE"],
#                                 'ed':row["ED"]}, ignore_index=True)
#             '''
#             fold_i += 1
#         else:
#             fold_i += 1
#
#         #Inititate base pair tabulation variables
#         bond_order = []
#         bond_count = 0
#
#         #Iterate through sequence to assign nucleotides to structure type
#         m = 0
#         for pair in structure:
#             if pair == '(':
#                 bond_count += 1
#                 bond_order.append(bond_count)
#                 m += 1
#             elif pair == ')':
#                 bond_order.append(bond_count)
#                 bond_count -= 1
#                 m += 1
#             elif pair == '.':
#                 bond_order.append(0)
#                 m += 1
#             else:
#                 print("Error1")
#
#         #Initiate base_pair list
#         base_pairs = []
#
#         #Create empty variable named test
#         test = ""
#
#         #Iterate through bond order
#         j = 0
#         while j < len(bond_order):
#         if bond_order[j] != 0:
#             test = bond_order[j]
#             base_pairs.append(j+1)
#             bond_order[j] = 0
#             j += 1
#             k = 0
#             while k < len(bond_order):
#                 if bond_order[k] == test:
#                     base_pairs.append(k+1)
#                     bond_order[k] = 0
#                     k += 1
#                 else:
#                     k += 1
#         else:
#             j += 1
#
#         #Iterate through "base_pairs" "to define bps
#         l = 0
#         while l < len(base_pairs):
#         lbp = base_pairs[l]
#         rbp = base_pairs[l+1]
#
#         lb = str(sequence[int(lbp)-1])
#         rb = str(sequence[int(rbp)-1])
#
#         lbp_coord = int(int(lbp)+int(row["Start"])-1)
#         rbp_coord = int(int(rbp)+int(row["Start"])-1)
#
#         bp_dict = bp_dict.append({'nuc_i':lb, 'i_bp':lbp_coord,
#                             'nuc_j':rb, 'j_bp':rbp_coord,
#                     'z':row["Z-score"], 'mfe':row["NativeMFE"],
#                             'ed':row["ED"]}, ignore_index=True)
#         bp_dict = bp_dict.append({'nuc_i':rb, 'i_bp':rbp_coord,
#                             'nuc_j':lb, 'j_bp':lbp_coord,
#                     'z':row["Z-score"], 'mfe':row["NativeMFE"],
#                             'ed':row["ED"]}, ignore_index=True)
#         l += 2
