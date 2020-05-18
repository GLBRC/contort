
# **COnTORT** 

## Purpose:

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

COnTORT has been tested on Linux (Ubuntu and CentOS), MacOS (>10.13), and Windows 10.

## Input : 

 - The GenBank (.gbff) for your organism, downloaded from NCBI.
             
You are encouraged to create a new directory each time you run this script.
## **THIS WILL FAIL IF RUN ON VERSION OF PYTHON BELOW 3**

#### Usage:

1. You can install using pip and then run:

        pip install contort
        contort -f gds_result.txt -gb genbank.gbff

2. You can install using Anaconda (anaconda.org/):
       
        conda install -c kevinmyers contort
        contort -f gds_result.txt -gb genbank.gff

3. You can download the git repository and run the original scripts:

        git clone https://github.com/GLBRC/contort.git
        python3.6 contort.py -f gds_result.txt -gb genbank.gbff

#### Requirements:

 - python 3
 - Python modules argparse, Bio, ftplib, functools, GEOparse, os, pandas, re, shutil, sys, subprocess, sys, tkinter, time, urllib

## Output : 

COnTORT_organized_transcriptomic_data.txt is the primary output.

Primary output is a single text file containing the gene annotations on the first five columns
    and the mean centered expression data present in the remaining columns. The first row
    contains the headers references the specific GEO IDs for each experiment. 
    Note that annotations no present in the GenBank file are listed with "N/A".
    
    Locus_Tag    Old_Locus_Tag    Gene_Name    Gene_Synonyms    Product    GSM_ID_1    GSM_ID_2
    RSP_0002     N/A              spbB         N/A              H-NS       12.0         4.0

    There will be two files:  one where the results are mean-centered for each experiment and one where no mean-centering is performed.
    
Subdirectories are created for all the files and to organize the directory:
        
#### Directories created:
	- geo                     - the downloaded GDS, GSE, GSM and GPL directories
	- GEOannotate_results     - the results of running GEOparse to organize the annotation and data
	- GeneOrf_match_output    - the results of matching the GEOparse results to the gene IDs from the GFF
	- mean_centered_results   - the mean centered expression data for each experiment
	- FTP_files               - files used to download the data from GEO via FTP
	- log_files               - all log files from each step as well as other saved files

## Outline of steps & commands used in pipeline:

#### runFindGEOAddresses( GDSfile ):

	Parse GEOfile
    
    Opens and parses the GDSfile result txt file from the search term provided by the user.
    
    Creates new files:
        - GEO_FTP_Addresses.txt     - all GDS, GSE, GSM, GPL addresses in the file
        - GEO_FTP_directories.txt   - GDS, GSE, GSM, GPL directories that will be downloaded


#### run_is_ftp_dir( ftp_handle, name, guess_by_extension=True ):

    QC for FTP
    
    Determines if an item listed on the FTP server is a directory or not by 
    looking for a "." in the fourth position. If it has that, it is nearly 
    always a file and not a directory.

#### run_make_parent_dir( fpath ):

    Make directories to match the FTP
    Creates the directories in the local directory to match the FTP directories

#### run_download_FTP_file( ftp_handle, name, dest, overwrite ):

    Copy the FTP files to the local directory
    Copy FTP files into the respective directories on the local directory

#### run_mirror_ftp_dir( ftp_handle, name, overwrite, guess_by_extension ):

    Replicate the directories
    Replicates a direcotry from the FTP server onto the local drive recusively

#### download_ftp_tree( ftp_handle, path, destination, overwrite=False, guess_by_extension=True ):

    Performs the actions to download files from the NCBI FTP
    
    Perform the actions
    Downloads an entire directory tree from an ftp server to the local destination
    Will NOT overwrite files if present in the local directory

#### Default FTP Settings for NCBI GEO:

	server = 'ftp.ncbi.nlm.nih.gov'
	user = 'anonymous'
	password = ''
	destination = user input
	sources of files = from the runFindGEOAddresses module

#### GEOannotate( ):

    Organize the GEO series files
    
    This is a replacement for R and GEOquery. This will organize the GEO data available.
    
    This first finds and copies all soft.gz files for all GSE (series) data in the 
    GEO download from the previous steps. Using GEOparse, the metadata for each gene is 
    collected and concatentated with the normalized data present for each series.
    
    The new files are written for downstream steps and the copied soft.gz files are deleted.

#### gffMatch( GBFF ):

    Match the gene annotations from the GenBank file to the organized expression data
    
    Using dictionaries created from the GenBank file for the orgainsm, this
    script will search for matches to the gene annotation in the GenBank file and
    retain only those data with matches. This will make new files for each GEOquery
    output with columns representing gene annotations and then the log2 normalized data
    from the GEO series files. The data are then mean centered and joined together
    with the gene annotations annotation as the key. All blanks are retained for consistency.
    This file is written and can be used in Excel or R for further analysis.

#### cleanUp( cwd ):

    Clean up the directory
    
    Organize the files into folders for a cleaner directory
