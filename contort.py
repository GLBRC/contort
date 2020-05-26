#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""COnTORT (COmprehensive Transcriptomic ORganizational Tool)

Description
-----------

COnTROT (COmprehensive Transcriptomic ORganizational Tool) is a program that
will download and organize all expression data in GEO related to a search result.
This will identify and download GEO GDS, GSE, GSM, and GPL
directories and files from NCBI FTP from a GDS search result
at NCBI. Then the series data will be organized to retain matches to any 
gene annotation present in the GenBank file provided. The data will be
organized, mean centered, and then joined using the gene annotation into a
single text file that can be easily manipulated or opened in Excel.

The results will be two files: one with the data mean centered and the other with no mean centering.

The user must provide the GenBank file to use for annotation from NCBI.

COnTORT will request the search term and use that to search NCBI GEO and download the NCBI GEO results.

Note that COnTORT can only download and combine normalized and processed data imported into NCBI GEO.

Notes
-----

    Must be run on scarcity. Requires Python3

    inputs : The GenBank (.gbff) for your organism, downloaded from NCBI.
    
    output : Primary output is a single text file containing the gene IDs on the first column
    and the organized expression data present in the remaining columns. The first row
    contains the headers references the specific GEO IDs for each experiment. 
    Note that annotations no present in the GenBank file are listed with "N/A".
    There will be two files:  one where the results are mean-centered for each experiment and one where no mean-centering is performed.
    
        Locus_Tag    Old_Locus_Tag    Gene_Name    Gene_Synonyms    Product    GSM_ID_1    GSM_ID_2
        RSP_0002     N/A              spbB         N/A              H-NS       12.0         4.0
    
    Subdirectories are created for all the files and to organize the directory:
        
    Directories created:
        - geo                     - the downloaded GDS, GSE, GSM and GPL directories
        - GEOannotate_results     - the results of running GEOparse to organize the annotation and data
        - GeneOrf_match_output    - the results of matching the GEOparse results to the gene IDs from the GFF
        - mean_centered_results   - the mean centered expression data for each experiment
        - FTP_files               - files used to download the data from GEO via FTP
        - log_files               - all log files from each step as well as other saved files
        
    method :
        Step 1) The NCBI GEO search result from the term provided by the user is downloaded to a GDS file and is parsed into two new files
    
        Step 2) The GEO data from the GDSfile is downloaded to the directory via FTP
    
        Step 3) GEOParse is run to organize the GEO data for which normalized data is available in NCBI GEO
    
        Step 4) The organized GEO Series file is matched to
                the gene annotation information and the resulting data are or are not mean
                centered within each experiment and combined into the final output files.
        
        Step 5) Clean up and orgainze the directory
    
    dependencies : 
        - python 3
        - Python modules argparse, Bio, ftplib, functools, GEOparse, os, pandas, re, shutil, sys, subprocess, sys, tkinter, time, urllib

    usage [standard]:
        python CONTORT.py
    
    usage [from pip install]:
        contort
                
"""
# Import required Python modules
import argparse
from Bio import SeqIO
from Bio import Entrez
import ftplib
import GEOparse
import os
from os import path
import pandas as pd
from pathlib import Path
import re
import shutil
import sys
import tkinter as tk
from tkinter import messagebox, Label, Button, Entry, W, Tk, Toplevel, filedialog, ttk
import time
import urllib.request as request


def runFindGEOAddresses( GDSfile ):
    """
    Parse GEOfile
    
    Opens and parses the GDSfile result txt file from the result of the search term
    provided by the user in the program.
    
    Creates new files:
        GEO_FTP_Addresses.txt - all GDS, GSE, GSM, GPL addresses in the file
        GEO_FTP_directories.txt - GDS, GSE, GSM, GPL directories that will be downloaded
    """
    # Parse the GDSfile
    file = open(GDSfile, 'r')
    filetext = file.read()
    file.close()
    gds_matches = re.findall("ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS\S*", filetext)
    with open('GEO_FTP_Addresses.txt', 'w') as f:
        for row in gds_matches:
            f.write("%s\n" % str(row))
    gds_number = len(gds_matches)
    gse_matches = re.findall("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE\S*", filetext)        
    with open('GEO_FTP_Addresses.txt', 'a+') as f:
        for row in gse_matches:
            f.write("%s\n" % str(row))
    gse_number = len(gse_matches)
    gsm_matches = re.findall("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM\S*", filetext)        
    with open('GEO_FTP_Addresses.txt', 'a+') as f:
        for row in gsm_matches:
            f.write("%s\n" % str(row))
    gsm_number = len(gsm_matches)
    gpl_matches = re.findall("ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL\S*", filetext)        
    with open('GEO_FTP_Addresses.txt', 'a+') as f:
        for row in gpl_matches:
            f.write("%s\n" % str(row))
    gpl_number = len(gpl_matches)
    
    print(f"Parsing the NCBI GEO search results file...\n")
    print(f"There are {gds_number} GDS (dataset) addresses in the file.")
    print(f"There are {gse_number} GSE (series) addresses in the file.")
    print(f"There are {gsm_number} GSM (samples) addresses in the file.")
    print(f"There are {gpl_number} GPL (platforms) addresses in the file.\n")
    
    print(f"Downloading the data from NCBI GEO now...\n")
    # Add addresses and directories to lists and write to files    
    addresses = []
    
    with open('GEO_FTP_Addresses.txt', 'r') as f:
        for line in f:
            line = line.rstrip()
            addresses.append(line)
                
    directories = []
    
    for line in addresses:
        line = line.replace('ftp://ftp.ncbi.nlm.nih.gov/', '/')
        directories.append(line)
        
    with open('GEO_FTP_Directories.txt', 'w') as f:
        for row in directories:
            f.write("%s\n" % str(row))
        
def run_is_ftp_dir(ftp_handle, name, guess_by_extension=True):
    """
    QC for FTP
    
    Determines if an item listed on the FTP server is a directory or not by 
    looking for a "." in the fourth position. If it has that, it is nearly 
    always a file and not a directory.
    """
    
    if guess_by_extension is True:
        if len(name) >= 4:
            if name[-4] == '.':
                return False            # check if the file is directory or file

    original_cwd = ftp_handle.pwd()     # remember the current working directory
    try:
        ftp_handle.cwd(name)            # try to set directory to new name
        ftp_handle.cwd(original_cwd)    # set it back to what it was
        return True
    
    except ftplib.error_perm as e:      # error warning
        return False
    
    except:
        return False


def run_make_parent_dir(fpath):
    """
    Make directories to match the FTP
    Creates the directories in the local directory to match the FTP directories
    """
    
    dirname = os.path.dirname(fpath)
    while not os.path.exists(dirname):    # if the directory does not exist on the local drive, create it
        try:
            os.makedirs(dirname)
            with open('FTP_download_log.txt','a') as log:
                log.write("created {0}".format(dirname)) # write to log file
                log.close()
        except:
            run_make_parent_dir(dirname)


def run_download_ftp_file(ftp_handle, name, dest, overwrite):
    """
    Copy the FTP files to the local directory
    Copy FTP files into the respective directories on the local directory
    """
    
    run_make_parent_dir(dest.lstrip("/"))
    if not os.path.exists(dest) or overwrite is True:  # copy the file from FTP to the same directory on the local drive
        try:
            with open(dest, 'wb') as f:
                ftp_handle.retrbinary("RETR {0}".format(name), f.write) # copy file to local directory
            with open('FTP_download_log.txt','a') as log:
                log.write("downloaded: {0}\n".format(dest)) # write log file
                log.close()
        except FileNotFoundError:                      # error warning
            with open('FTP_download_log.txt','a') as log:
                log.write("FAILED: {0}".format(dest)) # write error in log file
                log.close()
    else:
        print("already exists: {0}".format(dest))


def run_mirror_ftp_dir(ftp_handle, name, overwrite, guess_by_extension):
    """
    Replicate the directories
    Replicates a direcotry from the FTP server onto the local drive recusively
    """
    
    for item in ftp_handle.nlst(name):
        if run_is_ftp_dir(ftp_handle, item, guess_by_extension):
            run_mirror_ftp_dir(ftp_handle, item, overwrite, guess_by_extension)
        else:
            run_download_ftp_file(ftp_handle, item, item, overwrite)


def download_ftp_tree(ftp_handle, FTPpath, destination, overwrite=False, guess_by_extension=True):
    """
    Performs the actions to download files from the NCBI FTP
    
    Perform the actions
    Downloads an entire directory tree from an ftp server to the local destination
    Will NOT overwrite files if present in the local directory
    """
   
    path = FTPpath.lstrip("/")
    original_directory = os.getcwd()    # remember working directory before function is executed
    os.chdir(destination)               # change working directory to ftp mirror directory
    run_mirror_ftp_dir(ftp_handle, path, overwrite, guess_by_extension)
    os.chdir(original_directory)        # reset working directory to what it was before function exec

def GEOannotate():
    """
    This will organize the GEO data available using GEOparse.
    
    This first finds and copies all soft.gz files for all GSE (series) data in the 
    GEO download from the previous steps. Using GEOparse, the metadata for each gene is 
    collected and concatentated with the normalized data present for each series.
    
    The new files are written for downstream steps and the copied soft.gz files are deleted.
    """
    for file in Path("./geo").rglob('*.soft.gz'):
        shutil.copy(file, './')
        
    to_remove = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("full.soft.gz") ]
    for each in to_remove:
        os.remove(each)
    
    to_remove2 = [ fn for fn in os.listdir(os.getcwd()) if fn.startswith("GDS") ]
    for each in to_remove2:
        os.remove(each)
    
    soft_files = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("soft.gz") ] 
    
    for each in soft_files: #need to work on skipping "each" if there is a keyError
        try:
            gse = GEOparse.get_GEO(filepath=each)
            samples = gse.pivot_samples("VALUE")
            samples.reset_index(level = 0, inplace = True)
            samples.rename(columns={'ID_REF':"ID"}, inplace = True)        
            
            for key,value in gse.gpls.items():
                GPL = key
            annotation = gse.gpls[GPL].table
            annotation.reset_index(level = 0, inplace = True)
    
            result = pd.merge(annotation, samples)
            gseID = each.split("_")[0]
            output = gseID + "_GEOannotate_results.txt"
            result.to_csv(output, sep = "\t", index = False)
        except:
            pass

    for each in soft_files:
        os.remove(each)

def gffMatch( GBFF ):
    """
    Match the gene annotations from the GenBank file to the organized expression data
    
    Using dictionaries created from the GenBank file for the orgainsm, this
    script will search for matches to the gene annotation in the GenBank file and
    retain only those data with matches. This will make new files for each GEOannotate
    output with columns representing gene annotations and then the log2 normalized data
    from the GEO series files. The data are then mean centered and joined together
    with the gene annotations annotation as the key. All blanks are retained for consistency.
    This file is written and can be used in Excel or R for further analysis.
    
    """
    combined = []
    combined2 = []
    
    # Make the gene annotation dictionary from the GenBank file
    
    for record in SeqIO.parse(GBFF, "genbank"):
        for f in record.features:
            if f.type == "CDS":
                if "locus_tag" in f.qualifiers:
                    locus_tag = f.qualifiers['locus_tag']
                else:
                    locus_tag = "N/A"
                if "gene" in f.qualifiers:
                    gene = f.qualifiers['gene']
                else:
                    gene = "N/A"
                if "old_locus_tag" in f.qualifiers:
                    old_locus_tag = f.qualifiers['old_locus_tag']
                else:
                    old_locus_tag = "N/A"
                if "product" in f.qualifiers:
                    product = f.qualifiers['product']
                else:
                    product = "N/A"
                if "gene_synonym" in f.qualifiers:
                    gene_synonym = f.qualifiers['gene_synonym']
                else:
                    gene_synonym = "N/A"
                    combined.append(f"{locus_tag}\t{old_locus_tag}\t{gene}\t{gene_synonym}\t{product}\n")
    
    for line in combined:
        line2 = re.sub("\]", "", line)
        line2 = re.sub("\[", "", line2)
        line2 = re.sub("['']", "", line2)
        combined2.append(line2)

    # Write gene annotations to file
    
    with open("gene_annotations.txt", 'w') as w:
        w.write("Locus_Tag\tOld_Locus_Tag\tGene_Name\tGene_Synonyms\tProduct\n")

    for each in combined2:
        with open("gene_annotations.txt", 'a') as w:
            w.write(each)
    
    combo = {}
    
    with open('gene_annotations.txt', 'r') as f:
        next(f)
        for line in f:
            locusTagID = line.split('\t')[0]
            oldLocusTagID = line.split('\t')[1]
            geneID = line.split('\t')[2]
            synonyms = line.split('\t')[3]
            products = line.split('\t')[4]
            productID = products.strip('\n')
            all_annotations = locusTagID + "\t" + oldLocusTagID + "\t" + geneID + "\t" + synonyms + "\t" + productID
            if locusTagID == "N/A":
                pass
            if locusTagID in combo:
                pass        
            else:
                combo[locusTagID] = []
                combo[locusTagID].append(all_annotations)
            if oldLocusTagID == "N/A":
                pass
            if oldLocusTagID in combo:
                pass
            else:
                combo[oldLocusTagID] = []
                combo[oldLocusTagID].append(all_annotations)
            if geneID == "N/A":
                pass
            if geneID in combo:
                pass
            else:
                combo[geneID] = []
                combo[geneID].append(all_annotations)
            if synonyms == "N/A":
                pass
            else:
                words = synonyms.split('; ')
                for each in words:
                    if each in combo:
                        pass
                    else:
                        combo[each] = []
                        combo[each].append(all_annotations)
    
    GEOannotate_results = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("GEOannotate_results.txt") ]

    # Correlate GEOannotate results with gene annotation dictionary
    
    for each in GEOannotate_results:
        with open(each, 'r', encoding = "utf-8", errors = 'ignore') as file:
            geoName = each.split('_')[0]
            firstLines = file.readlines()
            namesOnly = firstLines.pop(0)
            if namesOnly.startswith("GSM"):
                pass
            if "GSM" not in namesOnly:
                pass
            else:
                names = namesOnly.rstrip().split('\t')
                indices = [i for i, elem in enumerate(names) if 'GSM' in elem]
                dataNames = namesOnly.rstrip().split('\t')[indices[0]:]
                dataNamesHeader = '\t'.join(dataNames)
                with open(geoName + '_output_results.txt', 'w') as f:
                    f.write(f'Experiment_ORF\t{dataNamesHeader}\n')
                    for each in firstLines:
                        annotation = each.rstrip().split('\t')[:indices[0]]
                        dataLines = each.rstrip().split('\t')[indices[0]:]
                        for k in annotation:
                                if k in combo:
                                    orfs = k
                                    data = dataLines
                                    dataOut = '\t'.join(data)
                                    f.write(f"{orfs}\t{dataOut}\n")
                                    break

                        
    GEO_output_results = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("output_results.txt") ]

    # Combine and write non-mean centered results file

    # Rename any results file that is less than 1 MB (no data)   

    for each in GEO_output_results:
        if os.path.getsize(each) < 1000:
            name = each.split('_')[0]
            os.rename(each, name+"_output_results_small.txt")
            
    GEO_output_results2 = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("output_results.txt") ]
    
    # Remove any duplicate values that appear in the mean center results files
    
    for f in GEO_output_results2:
        df_name = f.split("_")[0]
        df = pd.read_csv(f, sep = "\t")
        df.sort_values("Experiment_ORF", inplace = True)
        df2 = df.drop_duplicates(subset = "Experiment_ORF", keep = 'first')
        outname = df_name + "_output_results_noDuplicates.txt"
        df2.to_csv(outname, sep = "\t", header = True, index = False)
        
    GEO_output_results_noDups = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("output_results_noDuplicates.txt") ]
    
    # Combine all the results files into the final data file
    
    outputCombo = {}  
    outputCombo['Locus_Tag\tOld_Locus_Tag\tGene_Name\tGene_Synonyms\tProduct'] = []
    for each in GEO_output_results_noDups:    
        with open(each, 'r') as g:
            firstLine = g.readlines()
            header = firstLine[0]
            header2 = header.rstrip('\n').split('\t')[1:] 
            outputCombo['Locus_Tag\tOld_Locus_Tag\tGene_Name\tGene_Synonyms\tProduct'].append(header2)
            for each in firstLine:
                name2 = each.split('\t')[0]
                data = each.rstrip('\n').split('\t')[1:]
                if name2 in combo:
                    name = str(combo[name2])
                    if name not in outputCombo:
                        outputCombo[name]=[]
                        outputCombo[name].append(data)
                    else:
                        outputCombo[name].append(data)
                        
    combined_output_final = []
    
    for key,val in outputCombo.items():
        val2 = [val for sublist in outputCombo[key] for val in sublist]
        val3 = "\t".join(val2)    
        combined_output_final.append(f"{key}\t{val3}\n")
        combined_output_final[1:] = sorted(combined_output_final[1:])
    
    combined_output_final2 = []    
    for line in combined_output_final:
        line2 = re.sub("\]", "", line)
        line2 = re.sub("\[", "", line2)
        line2 = re.sub("['']", "", line2)
        line2 = line2.replace('\\t', '\t')
        combined_output_final2.append(line2)

    # Write COnTORT results file (Non-Mean Centered)
        
    final_original_output_name = "COnTORT_organized_non-mean_centered_transcriptomic_data.txt"
    
    with open(final_original_output_name, 'w') as f:
        for each in combined_output_final2:
            f.write(each)    
    
    # Mean center the organized results for each experiment to make comparisons easier
            
    GEO_output_results = []
    
    GEO_output_results = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("output_results.txt") ]
    
    for each in GEO_output_results:
        geoName = each.split('_')[0]
        df = pd.read_csv(each, sep = "\t", index_col = 0)
        df_MeanCenter = df.sub(df.mean(axis=1), axis=0)
        output_name = geoName + "_meanCenter_results.txt"
        df_MeanCenter.to_csv(output_name, sep = "\t")
        
    df_list_all = []
    
    df_list_all = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("meanCenter_results.txt") ]
    
    # Delete any mean centered file that is less than 1 MB (no data)   
    
    for each in df_list_all:
        if os.path.getsize(each) < 1000:
            os.remove(each)
    
    df_list_updated = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("meanCenter_results.txt") ]

    # Remove any duplicate values that appear in the mean center results files
        
    for f in df_list_updated:
        df_name = f.split("_")[0]
        df = pd.read_csv(f, sep = "\t")
        df.sort_values("Experiment_ORF", inplace = True)
        df2 = df.drop_duplicates(subset = "Experiment_ORF", keep = 'first')
        outname = df_name + "_meanCenter_results_noDuplicates.txt"
        df2.to_csv(outname, sep = "\t", header = True, index = False)

    
    df_list_final = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("meanCenter_results_noDuplicates.txt") ]

    # Combine all the mean centered results files into the final data file

    combo2 = {}  
    combo2['Locus_Tag\tOld_Locus_Tag\tGene_Name\tGene_Synonyms\tProduct'] = []
    for each in df_list_final:    
        with open(each, 'r') as g:
            firstLine = g.readlines()
            header = firstLine[0]
            header2 = header.rstrip('\n').split('\t')[1:] 
            combo2['Locus_Tag\tOld_Locus_Tag\tGene_Name\tGene_Synonyms\tProduct'].append(header2)
            for each in firstLine:
                name2 = each.split('\t')[0]
                data = each.rstrip('\n').split('\t')[1:]
                if name2 in combo:
                    name = str(combo[name2])
                    if name not in combo2:
                        combo2[name]=[]
                        combo2[name].append(data)
                    else:
                        combo2[name].append(data)
                        
    combined_final = []
    
    for key,val in combo2.items():
        val2 = [val for sublist in combo2[key] for val in sublist]
        val3 = "\t".join(val2)    
        combined_final.append(f"{key}\t{val3}\n")
        combined_final[1:] = sorted(combined_final[1:])
    
    combined_final2 = []    
    for line in combined_final:
        line2 = re.sub("\]", "", line)
        line2 = re.sub("\[", "", line2)
        line2 = re.sub("['']", "", line2)
        line2 = line2.replace('\\t', '\t')
        combined_final2.append(line2)

    # Write COnTORT mean-centered results file
        
    final_output_name = "COnTORT_organized_mean_centered_transcriptomic_data.txt"
    
    with open(final_output_name, 'w') as f:
        for each in combined_final2:
            f.write(each)
            
def cleanUp( cwd ):
    """
    Clean up the dictory
    
    Organize the files into folders for a cleaner directory
    """
    cwd = os.getcwd() + "/"    
    # Organize GEOquery results files
    os.mkdir( "GEOannotate_results" )
    gffDir = cwd + "/GEOannotate_results/"
    [ os.rename( (cwd + fn), (gffDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("GEOannotate_results.txt") ]
    # Organize GFF match output result files
    os.mkdir( "GeneOrf_match_output" )
    inputContigsDir = cwd + "/GeneOrf_match_output/"
    [ os.rename( (cwd + fn), (inputContigsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("output_results.txt") ]
    [ os.rename( (cwd + fn), (inputContigsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("output_results_small.txt") ]
    [ os.rename( (cwd + fn), (inputContigsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("output_results_noDuplicates.txt") ]
    # Organize mean centered result files
    os.mkdir( "mean_centered_results" )
    patserResultsDir = cwd + "/mean_centered_results/"
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("meanCenter_results.txt") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("meanCenter_results_noDuplicates.txt") ]
    # Organize log files
    os.mkdir( "log_files" )
    patserResultsDir = cwd + "/log_files/"
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("log.txt") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("gene_annotations.txt") ]
    # Organize GEO FTP files
    os.mkdir( "FTP_files" )
    patserResultsDir = cwd + "/FTP_files/"
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.startswith("GEO_FTP") ]

def main():
    
    cmdparser = argparse.ArgumentParser(description="Find and Download NCBI files from GDS file downloaded from NCBI GEO, process the data, mean center. and combine the data by GenBank annotation.", usage='%(prog)s -f <gds_result.txt> -gb <genbank_file.gbff>', prog='CONTORT.py' )
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL',  help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    cwd = os.getcwd()
    
    start = time.time() #start timer to time how long script takes to run
    
    if cmdResults['DETAIL']:
        print("\nScript Name: COnTORT")
        print("\tThis script was designed my Kevin Myers (kmyers2@wisc.edu).")
        print("\nPurpose : Download NCBI GEO Data from NCBI FTP. The data are then organized, combined with annotation from supplied GenBank files, mean centered, and concatenated to a single file.")
        print("\nInput : Using the GUI window, enter the search term for NCBI GEO and select the GenBank annotation file to use.")
        print("\nThis requires the the following Python modules: argparse, Bio, ftplib, functools, GEOparse, os, pandas, re, shutil, sys, subprocess, sys, tkinter, time, urllib.")
        print("\tDEFAULT Server:  ftp.ncbi.nlm.nih.gov")
        print("\tDEFAULT User:  anonymous")
        print("\tDEFAULT password:  \n")
        print("To run enter: contort\n")
        print("Then use the window that appears to enter the required inforamtion\n")
        print("See Kevin Myers (kmyers2@wisc.edu) for problems with COnTORT.\n\n")
        sys.exit(1)
        
    print("\nThanks for using COnTORT! Let's get started...\n")
    
    #TKinter class to make it easier for the user to enter the NCBI GEO search term and select the GenBank file to use
    
    class GetInfo:
        def __init__(self, master):
            self.master = master
            master.title("COnTORT")
            self.GenBank_File = None
            self.file_only = None
            
            self.label = Label(master, text = "Please enter the information for COnTORT to run.", font=("Arial Bold", 15))
            self.label.grid(column=0, row=0)
            
            self.GEOsearch_label = Label(master, text = "GEO Search Term:", font=("Arial Bold", 13)).grid(row = 1, column = 0, sticky = W)
            self.GEOsearch_explanation = Label(master, text = "Enter the phrase above you would like to search NCBI GEO for experiments in").grid(row = 3, column = 0, sticky = W)
            self.GEOsearch_explanation2 = Label(master, text = "High-throughput RNA Sequencing, Expression Microarrays, and Tiled Expression Microarrays.").grid(row = 4, column = 0, sticky = W)
            self.GEOsearch_explanation3 = Label(master, text = "COnTORT will analyze all available data on NCBI GEO matching the search results.").grid(row = 5, column = 0, sticky = W)
            self.GEOsearch_explanation3 = Label(master, text = "The more specific search term, the fewer results will be obtained.").grid(row = 6, column = 0, sticky = W)
            self.GEOsearch_explanation4 = Label(master, text = "  ").grid(row = 7, column = 0, sticky = W)
            self.GEOsearch_explanation5 = Label(master, text = "  ").grid(row = 8, column = 0, sticky = W)

            self.GenBankFileLabel = Label(master, text = "Open the GenBank File to Use", font=("Arial Bold", 13)).grid(row = 9, column = 0, sticky = W)
            self.GEO_term = Entry(master, width = 48)
            openGenBankFileCommand = master.register(self.GenBankFileOpen)
            self.genbank_button = ttk.Button(master, text = "Select GenBank File",
                                         command = openGenBankFileCommand).grid(row=10, column = 0, sticky = W, pady=10)
            self.GenBank_explanation = Label(master, text = "Use the button above to open the GenBank file (.gbff) for the annotation you would like to use in COnTORT.").grid(row = 11, column = 0, sticky = W)
            self.GenBank_explanation2 = Label(master, text = "GenBank annotation files can be downloaded for the specific strain and species").grid(row = 12, column = 0, sticky = W)
            self.GenBank_explanation3 = Label(master, text = "from https://www.ncbi.nlm.nih.gov/genome/.").grid(row = 13, column = 0, sticky = W)
            self.GenBank_explanation3 = Label(master, text = "Use the GenBank annotation file to match the NCBI GEO search of interest.").grid(row = 14, column = 0, sticky = W)
            self.GEOsearch_explanation4 = Label(master, text = "  ").grid(row = 15, column = 0, sticky = W)
            self.GEOsearch_explanation5 = Label(master, text = "  ").grid(row = 16, column = 0, sticky = W)
            
            self.GEO_term.grid(row = 1, column = 0)
            
            runContortCommand = master.register(self.getInput)
            exitContortCommand = master.register(self.getCancel)
            
            self.run_button=ttk.Button(master, text = "Run",
                       command = runContortCommand).grid(row = 20, column = 0, sticky = W)
            
            self.cancel_button=ttk.Button(master, text = "Cancel",
                       command = exitContortCommand).grid(row = 20, column = 1, sticky = W)
            
        def GenBankFileOpen(self):
            filename = filedialog.askopenfilename(title = "Open the GenBank file", filetypes = (("GenBank Files", "*.gbff"),("All Files", "*")))
            if filename: 
                self.GenBank_File = filename
                self.file_only = self.GenBank_File.split('/')[-1]
            
        def getInput(self):
            self.GEO_Search_Term = self.GEO_term.get()
            #self.GenBank_File_use = self.GenBank_File
            #self.GenBank_File_use
            if not self.GenBank_File:
                self.warning_window = tk.messagebox.showerror('Error', 'Please select a GenBank file to use.')
            else:
                self.close_box_window = tk.messagebox.askokcancel('Running COnTORT', "COnTORT will search and download the results for the following search term:\n{}\n\nCOnTORT will use the GenBank file: {}.".format(self.GEO_term.get(), self.file_only), default = 'ok')
                if self.close_box_window == True:
                    root.destroy()
                else:
                    return
        
        def getCancel(self):
            
            self.MsgBox_window = tk.messagebox.askokcancel("Exit COnTORT", "Are you sure you want to exit COnTORT?", icon = "warning", default = 'cancel')
            if self.MsgBox_window == False:
                return
            else:
                root.destroy()

    #Use TKinter to open a window to enter the NCBI GEO Search Term and select the GenBank File to use.

    root = Tk()
    root.geometry('800x500')
    contort_gui = GetInfo(root)

    root.mainloop()
    
    GDSsearch = contort_gui.GEO_Search_Term
    try:
        GBFF = contort_gui.file_only
    except:
        exit
    
    #TO GENERATE THE GDS RESULT FILE USING ENTREZ from BIOPYTHON:

    Entrez.email='kmyers2@wisc.edu'
    Entrez.api_key = "c47154faf0e8e818a313b184dc893b7b5a08"
    idlist = []
    handle = Entrez.esearch(db="gds", term=GDSsearch, retmax = '100000')
    result = Entrez.read(handle)
    for each in result['IdList']:
        idlist.append(each)
    for each in idlist:
        test = Entrez.esummary(db="gds", id=each)
        record = Entrez.read(test)
        for entry in record:
            if "Expression profiling by array" in entry['gdsType']:
                with open('GDS_result.txt', 'a') as f:
                    f.write(f"{entry['title']}\n{entry['summary']}\nOrganism:\t{entry['taxon']}\nType:\t{entry['gdsType']}\nPlatform:\tGPL{entry['GPL']}\t\t{entry['n_samples']} Samples\nFTP download:\t{entry['FTPLink']}\nSeries\t\tAccession:\t{entry['Accession']}\tID:\t{entry['Id']}\n\n")
            elif "Expression profiling by high throughput sequencing" in entry['gdsType']:
                with open('GDS_result.txt', 'a') as f:
                    f.write(f"{entry['title']}\n{entry['summary']}\nOrganism:\t{entry['taxon']}\nType:\t{entry['gdsType']}\nPlatform:\tGPL{entry['GPL']}\t\t{entry['n_samples']} Samples\nFTP download:\t{entry['FTPLink']}\nSeries\t\tAccession:\t{entry['Accession']}\tID:\t{entry['Id']}\n\n")
            elif "Expression profiling by genome tiling array" in entry['gdsType']:
                with open('GDS_result.txt', 'a') as f:
                    f.write(f"{entry['title']}\n{entry['summary']}\nOrganism:\t{entry['taxon']}\nType:\t{entry['gdsType']}\nPlatform:\tGPL{entry['GPL']}\t\t{entry['n_samples']} Samples\nFTP download:\t{entry['FTPLink']}\nSeries\t\tAccession:\t{entry['Accession']}\tID:\t{entry['Id']}\n\n")    

    GDSfile = 'GDS_result.txt'
    
    runFindGEOAddresses( GDSfile )
        
    server = 'ftp.ncbi.nlm.nih.gov'
    user = 'anonymous'
    password = ''
    destination = '.'
    ftp = ftplib.FTP(server, user, password)
      
    with open('GEO_FTP_Directories.txt', 'r') as f:
        for line in f:
            FTPpath = line.rstrip()
            download_ftp_tree(ftp, FTPpath, destination)

    print("Organizing and combining the transcriptomic data downloaded from NCBI GEO...\n") 
    
    GBFF = contort_gui.file_only
    
    GEOannotate()
    gffMatch( GBFF )
    cleanUp( cwd )
    
    # end timer and do math and report how long the script took to run
    end = time.time()
    total_time = round(end - start, 2)
    total_time_min = round(total_time/60, 2)
    total_time_hours = round(total_time/60/60, 2)
    print(f"\n\nIt took {total_time_hours} hours ({total_time_min} minutes) to download, organize, and combine the data.\n")
    print(f"We hope you enjoyed using COnTORT! Please e-mail Kevin Myers (kmyers2@wisc.edu) if you have any questions.\n")

if __name__ == "__main__":
    main()