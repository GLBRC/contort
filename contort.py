#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""COnTORT (COmprehensive Transcriptomic ORganizational Tool)

Description
-----------

COnTROT (COmprehensive Transcriptomic ORganizational Tool) is a program that
will download and organize all expression data in GEO related to a search result,
commonly an organism. This will identify and download GEO GDS, GSE, GSM, and GPL
directories and files from NCBI FTP when provided with a downloaded GDS result
from NCBI. Then the series data will be organized to retain matches to any 
gene annotation present in the GenBank file provided. The data will be
organized, mean centered, and then joined using the gene annotation into a
single text file that can be easily manipulated or opened in Excel.

To generate the initial file, perform a search in NCIB for GEO DataSets. Then in the
upper right hand corner of the search results, select 'Send to:' and choose
'File' and then click 'Create File'. This will contain a summary of your search
results. This script will search this file for series indicators and download
the respective files and directories. You will also need the GenBank file for your organism.

Notes
-----

    Must be run on scarcity. Requires Python3

    inputs : A gds_result text file downloaded from NCBI.
             The GenBank (.gbff) for your organism, downloaded from NCBI.
    
    output : Primary output is a single text file containing the gene IDs on the first column
    and the mean centered expression data present in the remaining columns. The first row
    contains the headers references the specific GEO IDs for each experiment. 
    Note that annotations no present in the GenBank file are listed with "N/A".
    
        Locus_Tag    Old_Locus_Tag    Gene_Name    Gene_Synonyms    Product    GSM_ID_1    GSM_ID_2
        RSP_0002     N/A              spbB         N/A              H-NS       12.0         4.0
    
    Subdirectories are created for all the files and to organize the directory:
        
    Directories created:
        - geo                     - the downloaded GDS, GSE, GSM and GPL directories
        - matrix_files            - the series_matrix.txt files downloaded from GEO
        - GEOquery_results        - the results of the GEOquery.r script
        - GFF_match_output        - the results of matching the GEOquery results to the gene IDs from the GFF
        - mean_centered_results   - the mean centered expression data for each experiment
        - FTP_files               - files used to download the data from GEO via FTP
        - log_files               - all log files from each step as well as other saved files
        
    method :
        Step 1) The user submitted GDSfile is parsed into two new files
    
        Step 2) The GEO data from the GDSfile is downloaded to the directory via FTP
    
        Step 3) GEOquery.r is run to organize the GEO Series data
    
        Step 4) The organized GEO Series file is matched to
                the gene annotation information and the resulting data are mean
                centered within each experiment and combined into the final output file.
        
        Step 5) Clean up and orgainze the directory
    
    dependencies : 
        - python 3
        - R libraries GEOquery and tidyverse
        - Python modles argparse, ftplib, functools, glob, gzip, io, os, pandas, re, shutil, subprocess, sys, time
        - GEO_annotate.r and GEOpatch.r Rscript in the same directory as COnTROT.py

    usage [standard]:
        python3.6 CONTORT.py -f gds_result.txt -gb genbank.gbff
                
        run with nohup and & to prevent timeouts
"""
# Import required Python modules
import argparse
from Bio import SeqIO
import ftplib
import glob
import gzip
import io
import os
import pandas as pd
import re
import shutil
import subprocess
import sys
import time
import urllib.request


def runFindGEOAddresses( GDSfile ):
    """
    Parse GEOfile
    
    Opens and parses the GDSfile result txt file provided by the user (downloaded
    from NCBI).
    
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


def download_ftp_tree(ftp_handle, path, destination, overwrite=False, guess_by_extension=True):
    """
    Performs the actions to download files from the NCBI FTP
    
    Perform the actions
    Downloads an entire directory tree from an ftp server to the local destination
    Will NOT overwrite files if present in the local directory
    """
   
    path = path.lstrip("/")
    original_directory = os.getcwd()    # remember working directory before function is executed
    os.chdir(destination)               # change working directory to ftp mirror directory
    run_mirror_ftp_dir(ftp_handle, path, overwrite, guess_by_extension)
    os.chdir(original_directory)        # reset working directory to what it was before function exec

def GEOquery():
    """
    Organize the GEO series files
    
    This first finds and copies all series_matrix.txt.gz files from the FTP download.
    These files are then unzipped and used to create the series.txt input file.
    Run the Rscript GEO_annotate.r - requires GEOpatch.r
    Runs GEOquery in R to download and combine the annotation with the log2 normalized
    data in the GSE files.
    This is then written to a new file for use in downstream steps.
    """
    for file in glob.glob('./**/*matrix.txt*', recursive = True):
        shutil.copy(file, '.')
        
    zipped = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("gz") ] # collect all zipped FASTQ files
    for each in zipped:
        out = each.split('.gz')[0]
        with gzip.open(each, 'r') as f_in, open(out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out) # unzip .gz files
            os.remove(each)

    series = []
    
    # Open directories file and get GSE entires. Write to list
    
    with open('GEO_FTP_Directories.txt', 'r') as f:
        for line in f:
            series_name = line.rstrip().split('/')[4]
            if series_name.startswith("GDS"):
                continue
            if series_name.startswith("GSE"):
                series.append(series_name)
            
    with open('series.txt', 'a+') as f:
        for row in series:
            f.write("%s\n" % str(row)) #save series entires to file
    
    # Run GEO_annotate.r to pull annotation files for each series entiry
    this_dir, this_filename = os.path.split(__file__)
    RSCRIPT_PATH = os.path.join(this_dir, "GEO_annotate.r")    
    cmd = [ 'Rscript', RSCRIPT_PATH ]
    output = subprocess.Popen(cmd, stderr=subprocess.PIPE).communicate()
    result = output[1].decode('utf-8')
    
    # write errors to log file
    
    with open('GEOquery_log.txt', 'a') as log:
        log.write(result)
        log.write("\n")

def gffMatch( GBFF ):
    """
    Match the gene annotations from the GenBank file to the organized expression data
    
    Using dictionaries created from the GenBank file for the orgainsm, this
    script will search for matches to the gene annotation in the GenBank file and
    retain only those data with matches. This will make new files for each GEOquery
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
    
    GEOquery_results_all = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("GEOquery_results.txt") ]

    # Remove any GEOquery results that are less than 1 MB (no data)
    
    for each in GEOquery_results_all:
        if os.path.getsize(each) < 900:
            os.remove(each)
    
    GEOquery_results = [ fn for fn in os.listdir(os.getcwd()) if fn.endswith("GEOquery_results.txt") ]

    # Correlate GEOquery results with gene annotation dictionary
    
    for each in GEOquery_results:
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

    # Mean center the organized results for each experiment to make comparisons easier
    
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

    # Remove any duplicate values that appear int he mean center results files
        
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

    # Write COnTORT results file
        
    final_output_name = "COnTORT_organized_transcriptomic_data.txt"
    
    with open(final_output_name, 'w') as f:
        for each in combined_final2:
            f.write(each)
            
def cleanUp( cwd ):
    """
    Clean up the dictory
    
    Organize the files into folders for a cleaner directory
    """
    cwd = os.getcwd() + "/"    
    # Organize downloaded series matrix files
    os.mkdir( "matrix_files" )
    originalContigDir = cwd + "/matrix_files/"
    [ os.rename( (cwd + fn), (originalContigDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("series_matrix.txt") ]
    [ os.rename( (cwd + fn), (originalContigDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("matrix.txt") ]
    # Organize GEOquery results files
    os.mkdir( "GEOquery_results" )
    gffDir = cwd + "/GEOquery_results/"
    [ os.rename( (cwd + fn), (gffDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("GEOquery_results.txt") ]
    # Organize GFF match output result files
    os.mkdir( "GeneOrf_match_output" )
    inputContigsDir = cwd + "/GeneOrf_match_output/"
    [ os.rename( (cwd + fn), (inputContigsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("output_results.txt") ]
    [ os.rename( (cwd + fn), (gffDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("RPKM.txt") ]
    # Organize mean centered result files
    os.mkdir( "mean_centered_results" )
    patserResultsDir = cwd + "/mean_centered_results/"
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("meanCenter_results.txt") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("meanCenter_results_noDuplicates.txt") ]
    # Organize log files
    os.mkdir( "log_files" )
    patserResultsDir = cwd + "/log_files/"
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("log.txt") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.startswith("series") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.startswith("gene_orf") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".bt2") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".fai") ]
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("sra_files.txt") ]
    # Organize GEO FTP files
    os.mkdir( "FTP_files" )
    patserResultsDir = cwd + "/FTP_files/"
    [ os.rename( (cwd + fn), (patserResultsDir + fn) ) for fn in os.listdir(cwd) if fn.startswith("GEO_FTP") ]

def main():
    
    cmdparser = argparse.ArgumentParser(description="Find and Download NCBI files from GDS file downloaded from NCBI GEO, process the data, mean center. and combine the data by GenBank annotation.", usage='%(prog)s -f <gds_result.txt> -gb <genbank_file.gbff>', prog='CONTORT.py' )
    cmdparser.add_argument('-f', '--file', action='store', dest='FILE', help = 'GDS_result.txt file downloaded from NCBI GEO.', metavar='')
    cmdparser.add_argument('-gb', '--genbank', action='store', dest='GBFF', help = 'Genbank file for the organism of interest.', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL',  help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    cwd = os.getcwd()
    
    start = time.time() #start timer to time how long script takes to run
    
    # if no args, print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
            
    if cmdResults['DETAIL']:
        print("\nScript Name: COnTORT.py")
        print("\tThis script was designed to run on Scarcity by kmyers on 2019-05-14.")
        print("\nPurpose : Bulk download of NCBI GEO Data from NCBI FTP. The data are then organized, combined with annotation from supplied GenBank files, mean centered, and concatenated.")
        print("\nInput : A text file containing the GDS result file downloaded from NCBI.")
        print("\nYou are encouranged to create a new directory in which to run this script.")
        print("\nThis requires the R libraries GEOquery and tidyverse to be installed.")
        print("Required Parameters: -f gds_result.txt\n")
        print("Required Parameters: -gb genbank_file.gbff\n")
        print("\tDEFAULT Server:  ftp.ncbi.nlm.nih.gov")
        print("\tDEFAULT User:  anonymous")
        print("\tDEFAULT password:  \n")
        print("To run enter:  python3.6 COnTORT.py -f gds_result.txt -gb genbank_file.txt\n")
        print("This will identify and download GEO GDS, GSE, GSM, and GPL directories and files from NCBI FTP when provided with a downloaded GDS result from NCBI. Then the series data will be organized to retain matches to the GenBank annotation (gene and ORF). The data will be organized, mean centered, and then joined using the gene annotations as IDs. \nTo generate the initial file, perform a search in NCIB for GEO DataSets. Then in the upper right hand corner of the search results, select 'Send to:' and choose 'File' and then click 'Create File'. This will contain a summary of your search results. This script will search this file for series indicators and download the respective files and directories.\n")
        print("See Kevin Myers (kmyers2@wisc.edu) for problems with COnTORT.\n\n")
        sys.exit(1)
        
    if cmdResults['FILE'] is not None:
        GDSfile = cmdResults['FILE']
    
    if cmdResults['GBFF'] is not None:
        GBFF = cmdResults['GBFF']
        
    if urllib.request.urlopen('http://www.google.com').getcode() == 0:
        print("\nUh oh! It looks like you don't have access to the interet. CoNTORT requires internet access to function. Please try again later when you have access to the internet.\n")
    else:
        pass
    
    if urllib.request.urlopen('http://www.apple.com').getcode() == 0:
        print("\nUh oh! It looks like you don't have access to the interet. CoNTORT requires internet access to function. Please try again later when you have access to the internet.\n")
    else:
        pass
    
    if urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/geo').getcode() == 0:
        print("\nUh oh! It looks like NCBI GEO is down. COnTORT cannot work without access to NCBI GEO. Please try again later when the webiste is available.\n")
    else:
        pass
    
    if urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/sra').getcode() == 0:
        print("\nUh oh! It looks like NCBI SRA is down. COnTORT cannot work without access to NCBI SRA. Please try again later when the webiste is available.\n")
    else:
        pass
    
    print("\nThanks for using COnTORT! Let's get started...\n")
    
    runFindGEOAddresses( GDSfile )
        
    server = 'ftp.ncbi.nlm.nih.gov'
    user = 'anonymous'
    password = ''
    destination = '.'
    ftp = ftplib.FTP(server, user, password)
      
    with open('GEO_FTP_Directories.txt', 'r') as f:
        for line in f:
            path = line.rstrip()
            download_ftp_tree(ftp, path, destination)

    print("Organizing and combining the transcriptomic data downloaded from NCBI GEO...\n") 
    
    GEOquery()
    gffMatch( GBFF )
    cleanUp( cwd )
    
    # end timer and do math and report how long the script took to run
    end = time.time()
    total_time = round(end - start, 2)
    total_time_min = round(total_time/60, 2)
    total_time_hours = round(total_time/60/60, 2)
    print(f"It took {total_time_hours} hours ({total_time_min} minutes) to download, organize, and combine the data.\n")
    print(f"We hope you enjoyed using COnTORT! Please e-mail Kevin Myers (kmyers2@wisc.edu) if you have any questions.\n")
                           
if __name__ == "__main__":
    main()