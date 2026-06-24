#!/usr/bin/python3
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

Contact: Ethan Coppenbarger - ecoppen@iastate.edu

Usage:
$ python3.6 ScanFold-Fold2.0.py fasta_filename --tsv tsv_filename [options]

"""

import argparse
import logging
import numpy as np
import os
import re
import shutil
import sys
import time
from Bio import SeqIO

from ScanFoldFunctions import *
import lib.fold.fold as fold

# import RNAstructure


def main(args):
    start_time = time.time()
    filename = args.filename
    tsv = args.tsv
    filter_value = int(args.filter)
    competition = int(args.competition)
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
    runfast = args.fast

    
    # parse string output of glob.glob into list
    if filename[0] == '[':
        filename = eval(filename)
    if tsv[0] == '[':
        tsv = eval(tsv)
    # TODO: change this, get working properly w/ list stuff
    filename = os.path.join(cwd, filename[0])
    # This sets strand to forward (but can be set up as an argument if needed)
    strand = 1
    for tsv_file in tsv:
        print(tsv_file)
        print(isinstance(tsv_file, str))
        os.chdir(cwd)
        # Read all lines from ScanFold-Scan file (including header)
        try:
            tsv_split = tsv_file.split(".tsv")
            outname = tsv_split[0]
        except:
            logging.inro("Input name should have .tsv extension")
            outname = tsv_file

        # handle case of input file being an absolute path
        outname = os.path.basename(outname)
        curr_dir = cwd
        if len(tsv) > 1:
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
        print(f"cwd: {cwd}")
        print(f"curr_dir: {curr_dir}")
        if folder_name is not None:
            try:
                if len(tsv) > 1:
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
                if len(tsv) > 1:
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
        # read tsv
        try:
            tsv_file.encode('utf-8').decode('utf-8')
            print("string encode/decode success")
            print(tsv_file)
        except UnicodeDecodeError:
            raise ValueError("string not in UTF-8")
        py_cur_cwd = os.getcwd()
        print(f"python cwd (actual): {py_cur_cwd}")
        print(f"python cwd (variable): {cwd}")
        tsv_file_path = os.path.join(cwd, tsv_file)
        base_pair_matrix = fold.BasePairMatrix(tsv_file_path)
        #base_pair_matrix.toCSV(mat_name)
        all_bps = base_pair_matrix.get_all_pairs()
        step_size = base_pair_matrix.getStepSize()
        seq_len = base_pair_matrix.getSequenceLength()
        win_size = base_pair_matrix.getWinSize()
        # read fasta
        # TODO: modify to work w/ multiple fastas
        # for loop is for when this function is added
        cur_record = SeqIO.parse(filename, "fasta")
        for r in cur_record:
            full_fasta_sequence = str(r.seq).upper()
            full_fasta_header = '>' + str(r.description)
            break
        #full_fasta_sequence = base_pair_matrix.getSequence()
        # Reset file names #
        dbn_file_path = outname + ".AllDBN_refold.txt"
        dbn_file_path1 = outname + ".NoFilter"
        dbn_file_path2 = outname + ".Zavg_-1"
        dbn_file_path3 = outname + ".Zavg_-2"
        dbn_file_path4 = outname + ".AllDBN.txt"
        
        '''
        Note that in ScanFold3 log files aren't generated in the same way as ScanFold2.
        Different algorithms mean that it wouldn't make sense to use the same format.
        '''
        #list of fold.BasePair
        if competition == 1:
            alg_start = time.process_time()
            final_partners = fold.greedy_approximation(base_pair_matrix)
            alg_end = time.process_time()
            elapsed_time = str(round((time.time() - start_time), 2))+"s"
            logging.info("time to run approximate matching algorithm: "+elapsed_time)
        elif competition == 0:
            # choose all base pairs with some value scanned
            final_partners = all_bps
            alg_start = time.process_time()
            alg_end = time.process_time()
            elapsed_time = str(round((time.time() - start_time), 2))+"s"
            logging.info("time to extract scanned pairs: "+elapsed_time)
        else:
            raise ValueError("Competition value not properly set")
        minus2_partners = filter_base_pairs(final_partners, -2.0)
        minus1_partners = filter_base_pairs(final_partners, -1.0)
        structure_extract_file = outname + ".structure_extracts.txt"
        refold_structure_extract_file = outname + ".refolded_structure_extracts.txt"
        #Determine start and end coordinate values
        with open(tsv_file_path, 'r') as ifile:
            ifile.readline()    # skip header
            tsv_line = ifile.readline().strip()
            start_coordinate = int(tsv_line.split()[0]) # read 1st column of second line
        with open(tsv_file_path, 'rb') as ifile:
            # following snippet walks back through file until it hits a newline
            try:
                ifile.seek(-2, os.SEEK_END)
                while ifile.read(1) != b'\n':
                    ifile.seek(-2, os.SEEK_CUR)
            except OSError:
                # one line file
                ifile.seek(0)
            tsv_line = ifile.readline().decode()
            end_coordinate = int(tsv_line.split()[0])

        # Write CT files
        if competition == 1:
            logging.info("Trying to write CT files with -c option")
            elapsed_time = str(round((time.time() - start_time), 2)) + "s"
            logging.info("Elapsed time: " + elapsed_time)
            logging.info("Writing CT files")
            ct_line_unfiltered = make_ct_lines(final_partners, full_fasta_sequence, strand, start_coordinate, end_coordinate)
            lines_to_ct(ct_line_unfiltered, str(dbn_file_path1)+".ct", sys.float_info.max, name)
            lines_to_ct(ct_line_unfiltered, str(dbn_file_path2)+".ct", -1.0, name)
            lines_to_ct(ct_line_unfiltered, str(dbn_file_path3)+".ct", -2.0, name)
            lines_to_ct(ct_line_unfiltered, (outname+".Zavg_0.ct"), 0.0, name)
            
            #list_to_ct(final_partners, (outname+".Zavg_nofilter.ct"), float(10), strand, name, start_coordinate, end_coordinate)
            #list_to_ct(final_partners, (outname+".Zavg_-2.ct"), float(-2), strand, name, start_coordinate, end_coordinate)
            #list_to_ct(final_partners, (outname+".Zavg_-1.ct"), float(-1), strand, name, start_coordinate, end_coordinate)
            if filter_value != -1:
                user_filter_partners = filter_base_pairs(final_partners, filter_value)
                print(filter_value)
                user_filter_dbn = outname + ".Zavg_" + str(filter_value)
                print(user_filter_dbn)
                #lines_to_ct(final_partners, str(user_filter_dbn)+".ct", filter_value, strand, name, start_coordinate)
                lines_to_ct(ct_line_unfiltered, str(user_filter_dbn)+".ct", filter_value, name)
                #makedbn_with_fasta(user_filter_dbn, filename, "Zavg"+str(filter_value))
                fold.make_dbn(user_filter_partners, str(full_fasta_sequence), str(full_fasta_header), str(user_filter_dbn+".dbn"))
                
            #makedbn_with_fasta(dbn_file_path1, filename, "NoFilter")
            #makedbn_with_fasta(dbn_file_path2, filename, "Zavg_-1")
            #makedbn_with_fasta(dbn_file_path3, filename, "Zavg_-2")
            print("len of partner lists")
            print(str(len(final_partners)))
            print(str(len(minus1_partners)))
            print(str(len(minus2_partners)))
            fold.make_dbn(final_partners, str(full_fasta_sequence), str(full_fasta_header), str(dbn_file_path1+".dbn"))
            fold.make_dbn(minus1_partners, str(full_fasta_sequence), full_fasta_header, (dbn_file_path2+".dbn"))
            fold.make_dbn(minus2_partners, str(full_fasta_sequence), full_fasta_header, (dbn_file_path3+".dbn"))
            igv_path_base = os.path.join(igv_path, outname)
            write_bp_from_list(final_partners, igv_path_base+".bp", start_coordinate, name)
            write_wig_list(final_partners, igv_path_base+".zavgs.wig", name, step_size, str("zscore"))
            write_wig_list(final_partners, igv_path_base+".mfe_avgs.wig", name, step_size, str("mfe"))
            write_wig_list(final_partners, igv_path_base+".ed_avgs.wig", name, step_size, str("ed"))
            write_bp_from_list(all_bps, igv_path_base+".ALL.bp", start_coordinate, name)
            
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
        # update past this point
        elif competition == 0:
            if args.bp_track is not None:
                #write_bp(best_bps, args.bp_track, start_coordinate, name, minz)
                write_bp_from_list(all_bps, args.bp_track, start_coordinate, name)
            else:
                write_bp_from_list(all_bps, igv_path_base+".ALL.bp", start_coordinate, name)
        else:
            raise ValueError("Competition value not properly set")
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
        # san checks
        length = len(full_fasta_sequence)
        length_st = len(full_filter_structure)
        if length != length_st:
            print(full_fasta_sequence)
            print(full_filter_structure)
            raise ValueError("Length of sequence and structure do not match")
        # extract structures
        extract_sequences = []  # list of strings
        extract_structures = [] # list of strings
        extract_coords = []     # list of tuples with start/end coords
        refold_structures = []
        motif_count = 0
        offset = 0  # set value to add unpaired nucleotides to 5' and 3' ends
        extracted_structure_list = []
        fold.structure_extract(full_fasta_sequence, 
                               full_filter_structure, 
                               extract_sequences, 
                               extract_structures, 
                               extract_coords)
        '''
        with open(structure_extract_file, 'w') as extract_file:
            for i, coord in enumerate(extract_coords):
                motif_count += 1
                motif_start_coord = coord[0] - offset
                motif_end_coord = coord[1] + offset
                if motif_start_coord < 0: 
                    raise ValueError(f"start coordinate ({motif_start_coord}) out of bounds for motif {motif_count}")
                if motif_end_coord > int(length): 
                    raise ValueError(f"end coordinate ({motif_end_coord}) out of bounds for motif {motif_count}")
                extracted_structure_list.append(ExtractedStructure(motif_count,
                                                                   extract_sequences[i],
                                                                   extract_structures[i],
                                                                   motif_start_coord,
                                                                   motif_end_coord))
                
                gff_attributes = f'motif_{motif_count};sequence={extract_sequences[i]};structure={extract_structures[i]}'
                extract_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, "ScanFold3", "RNA_sequence_secondary_structure", str(motif_start_coord), str(motif_end_coord), ".", ".", ".", gff_attributes))
        '''

        """ Obtain statistics for motifs """
        zscore_total = []
        pvalue_total = []
        mfe_total = []
        ed_total = []
        with open(refold_structure_extract_file, "w", newline='\n') as extract_file:
            for idx, sequence in enumerate(extract_sequences):
                motif_count += 1
                if not runfast:
                    fc = RNA.fold_compound(sequence)  # creates "Fold Compound" object
                    fc.hc_add_from_db(extract_structures[idx])
                    fc.pf()  # performs partition function calculations
                    frag_q = (RNA.pf_fold(sequence))  # calculate partition function "fold" of fragment
                    (MFE_structure, mfe_calc) = fc.mfe()  # calculate and define variables for mfe and structure
                    mfe_calc = round(mfe_calc, 2)
                    mfe_total.append(mfe_calc)
                    (centroid, distance) = fc.centroid()  # calculate and define variables for centroid
                    ed_calc = round(fc.mean_bp_distance(), 2)  # this calculates ED based on last calculated partition function
                    ed_total.append(ed_calc)
                    seqlist = []  # creates the list we will be filling with sequence fragments
                    seqlist.append(sequence)  # adds the native fragment to list
                    scrambled_sequences = scramble(sequence, 100, type)
                    seqlist.extend(scrambled_sequences)
                    energy_list = energies(seqlist, temperature, algo)
                    try:
                        zscore = round(zscore_function(energy_list, 100), 2)
                    except:
                        zscore = zscore_function(energy_list, 100)
                    zscore_total.append(zscore)
                    pvalue = round(pvalue_function(energy_list, 100), 2)
                    pvalue_total.append(pvalue)
                    # ps_title = f"motif_{motif_num} coordinates {es.i} - {es.j}"
                    es_dbn_path = f"{extract_path}/{name}_motif_{motif_count}.dbn"
                    with open(es_dbn_path, 'w') as es_dbn:
                        es_dbn.write(f">{name}_motif_{motif_count}_coordinates:{extract_coords[idx][0]}-{extract_coords[idx][1]}\n{sequence}\n{extract_structures[idx]}")
                    #dbn2ct(es_dbn_path)
                    es_dbn_outfile = es_dbn_path[:-3] + "ct"
                    fold.dbn_to_ct(es_dbn_path, es_dbn_outfile)
                    mot_struc_count_path = f"{name}_motif_{motif_count}.ct"
                    os.rename(os.path.join(extract_path, mot_struc_count_path),
                              os.path.join(inforna_path, mot_struc_count_path))

                    # Create postscript files
                    ps_fname = f"motif_{motif_count}.ps"
                    ps_loc = os.path.join(extract_path, ps_fname)
                    RNA.PS_rna_plot_a(sequence, MFE_structure, ps_loc, '',
                                      '')

                    # Set extracted structures up as GFF format
                    gff_attributes = f'motif_{motif_count};sequence={sequence};structure={extract_structures[idx]};refoldedMFE={MFE_structure};MFE(kcal/mol)={str(mfe_calc)};z-score={str(zscore)};ED={str(ed_calc)}'
                    extract_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, "ScanFold3", "RNA_sequence_secondary_structure", extract_coords[idx][0], extract_coords[idx][1], ".", ".", ".", gff_attributes))
                else:
                    es_dbn_path = f"{extract_path}/{name}_motif_{motif_count}.dbn"
                    with open(es_dbn_path, 'w') as es_dbn:
                        es_dbn.write(f">{name}_motif_{motif_count}_coordinates:{extract_coords[idx][0]}-{extract_coords[idx][1]}\n{sequence}\n{extract_structures[idx]}")
                    #dbn2ct(es_dbn_path)
                    es_dbn_outfile = es_dbn_path[:-3] + "ct"
                    fold.dbn_to_ct(es_dbn_path, es_dbn_outfile)
                    mot_struc_count_path = f"{name}_motif_{motif_count}.ct"
                    os.rename(os.path.join(extract_path, mot_struc_count_path),
                              os.path.join(inforna_path, mot_struc_count_path))

                    # Create postscript files
                    ps_fname = f"motif_{motif_count}.ps"
                    ps_loc = os.path.join(extract_path, ps_fname)
                    RNA.PS_rna_plot_a(sequence, extract_structures[idx], ps_loc, '',
                                      '')

                    # Set extracted structures up as GFF format
                    gff_attributes = f'motif_{motif_count};sequence={sequence};structure={extract_structures[idx]}'
                    extract_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, "ScanFold3", "RNA_sequence_secondary_structure", extract_coords[idx][0], extract_coords[idx][1], ".", ".", ".", gff_attributes))
            '''
            for motif in extracted_structure_list:
                if not runfast:
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
                    with open(es_dbn_path, 'w') as es_dbn:
                        es_dbn.write(f">{name}_motif_{motif.structure_count}_coordinates:{motif.i}-{motif.j}\n{motif.sequence}\n{motif.structure}")
                    dbn2ct(es_dbn_path)
                    mot_struc_count_path = f"{name}_motif_{motif.structure_count}.ct"
                    os.rename(os.path.join(extract_path, mot_struc_count_path),
                              os.path.join(inforna_path, mot_struc_count_path))

                    # Create postscript files
                    ps_fname = f"motif_{motif.structure_count}.ps"
                    ps_loc = os.path.join(extract_path, ps_fname)
                    RNA.PS_rna_plot_a(frag, MFE_structure, ps_loc, '',
                                      '')

                    # Set extracted structures up as GFF format
                    gff_attributes = f'motif_{motif.structure_count};sequence={motif.sequence};structure={str(motif.structure)};refoldedMFE={str(MFE_structure)};MFE(kcal/mol)={str(mfe_calc)};z-score={str(zscore)};ED={str(ed_calc)}'
                    extract_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), "ScanFold3", "RNA_sequence_secondary_structure", str(motif.i), str(motif.j), ".", ".", ".", gff_attributes))
                else:
                    es_dbn_path = f"{extract_path}/{name}_motif_{motif.structure_count}.dbn"
                    with open(es_dbn_path, 'w') as es_dbn:
                        es_dbn.write(f">{name}_motif_{motif.structure_count}_coordinates:{motif.i}-{motif.j}_zscore={zscore}\n{frag}\n{MFE_structure}")
                    dbn2ct(es_dbn_path)
                    mot_struc_count_path = f"{name}_motif_{motif.structure_count}.ct"
                    os.rename(os.path.join(extract_path, mot_struc_count_path),
                              os.path.join(inforna_path, mot_struc_count_path))

                    # Create postscript files
                    ps_fname = f"motif_{motif.structure_count}.ps"
                    ps_loc = os.path.join(extract_path, ps_fname)
                    RNA.PS_rna_plot_a(frag, MFE_structure, ps_loc, '',
                                      '')

                    # Set extracted structures up as GFF format
                    gff_attributes = f'motif_{motif.structure_count};sequence={motif.sequence};structure={str(motif.structure)};refoldedMFE={str(MFE_structure)};MFE(kcal/mol)={str(mfe_calc)};z-score={str(zscore)};ED={str(ed_calc)}'
                    extract_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), "ScanFold3", "RNA_sequence_secondary_structure", str(motif.i), str(motif.j), ".", ".", ".", gff_attributes))
            '''        
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

    # scan arguments
    parser.add_argument('filename',  nargs="+",
                        help='Input filename')
    parser.add_argument('-s', '--step', type=int, default=1,
                        help='Step size; default = 1')
    parser.add_argument('-w', '--window', type=int, default=120,
                        help='Window size; default = 120')
    parser.add_argument('-r', type=int, default=100,
            help='Number of randomizations for background shuffling; default = 100')
    parser.add_argument('--algorithm', type=str, default="rnafold",
            help='Folding algorithm used; rnafold, rnastructure, mxfold')

    # fold arguments
    parser.add_argument('--tsv',  nargs="+",
                        help='Input tsv name from ScanFold-Scan (ScanFoldFold.py only)')
    parser.add_argument('-f', '--filter', type=int, default=-1,
                        help='Filter value')
    parser.add_argument('-c', '--competition', type=int, default=1,
                        help='Competition (1 for disallow competition, 0 for allow; 1 by default)')
    parser.add_argument('--id', type=str, default = "UserInput",
                        help='Name or ID of sequence being analyzed. Default "UserInput"')
    parser.add_argument('--global_refold', action='store_true', default=False,
                        help='Global refold option. Refold full sequence using Zavg <-1 and <-2 base pairs')
    parser.add_argument('--folder_name',  type=str,
                        help='Name of output folder (defaults to header name or date/time)')
    parser.add_argument('--extract', type=int, default='1',
                        help='Extract structures from minus 1 or minus 2 dbn file (2 or 1); Default = 1')
    parser.add_argument('--es_path', type=str, default = "extracted_structures",
                        help='Name of extracted structures file')
    parser.add_argument('--igv_path', type=str, default = "igv_files",
                        help='Name of IGV file')
    parser.add_argument('--inforna_path', type=str, default = "inforna_structures",
                        help='Name of inforna file')
    parser.add_argument('--fast', action='store_true', default=False, 
                        help="skip over some steps in ScanFold-Fold, primarily refolding motifs, for faster runtime")
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
