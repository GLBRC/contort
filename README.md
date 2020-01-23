
# **COnTORT** 

## Purpose:

COnTROT (COmprehensive Transcriptomic ORganizational Tool) is a program that
will download and organize all expression data in GEO related to a search result,
commonly an organism. This will identify and download GEO GDS, GSE, GSM, and GPL
directories and files from NCBI FTP when provided with a downloaded GDS result
from NCBI. Then the series data will be organized to retain matches to the GenBank
annotation. The data will be organized, mean centered, and then
joined using the gene name and ORF as IDs into a single text file that can be
easily manipulated or opened in applications such as Excel.

To generate the initial file, perform a search in NCIB for GEO DataSets. Then in the
upper right hand corner of the search results, select 'Send to:' and choose
'File' and then click 'Create File'. This will contain a summary of your search
results. This script will search this file for series indicators and download
the respective files and directories. You will also need the GenBank file for your organism.

## Input : 

 - A gds_result text file downloaded from NCBI.

 - The GenBank (.gbff) for your organism, downloaded from NCBI.
             
You are encouraged to create a new directory each time you run this script.
## **THIS WILL FAIL IF RUN ON VERSION OF PYTHON BELOW 3**

#### Required Parameters:
	
	-f gds_result.txt
	-gb genebank_file.gbff

#### Usage:

	contort -f gds_result.txt -gb genbank.gbff

#### Requirements:

 - python 3
 - R libraries GEOquery and tidyverse
 - Python modles argparse, ftplib, functools, glob, gzip, io, os, pandas, re, shutil, subprocess, sys, time

To install the R libraries:

    install.packages("tidyverse")

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("GEOquery")

## Output : 

Primary output is a single text file containing the gene IDs on the first column
    and the mean centered expression data present in the remaining columns. The first row
    contains the headers references the specific GEO IDs for each experiment. 
    Note that annotations no present in the GenBank file are listed with "N/A".
    
    Locus_Tag    Old_Locus_Tag    Gene_Name    Gene_Synonyms    Product    GSM_ID_1    GSM_ID_2
    RSP_0002     N/A              spbB         N/A              H-NS       12.0         4.0
    
Subdirectories are created for all the files and to organize the directory:
        
#### Directories created:
	- geo                     - the downloaded GDS, GSE, GSM and GPL directories
	- matrix_files            - the series_matrix.txt files downloaded from GEO
	- GEOquery_results        - the results of the GEOquery.r script
	- GFF_match_output        - the results of matching the GEOquery results to the gene IDs from the GFF
	- mean_centered_results   - the mean centered expression data for each experiment
	- FTP_files               - files used to download the data from GEO via FTP
	- log_files               - all log files from each step as well as other saved files

## Outline of steps & commands used in pipeline:

#### runFindGEOAddresses( GDSfile ):

	Parse GEOfile
    
    Opens and parses the GDSfile result txt file provided by the user (downloaded
    from NCBI).
    
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

#### GEOquery( ):

	Organize the GEO series files
    
    This first finds and copies all series_matrix.txt.gz files from the FTP download.
    These files are then unzipped and used to create the series.txt input file.
    Run the Rscript GEO_annotate.r - requires GEOpatch.r
    Runs GEOquery in R to download and combine the annotation with the log2 normalized
    data in the GSE files.
    This is then written to a new file for use in downstream steps.

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
