# **COnTORT Information**

### Detailed Description of COnTORT 

COnTORT performs several steps to organize the information retrieved from NCBI GEO. The user defined NCBI GEO serch query results in a downloaded results file that is parsed and an additional file is produced containing this information. GEO_FTP_Addresses.txt contains all the FTP addresses for all GPL (GEO Platform), GSE (GEO Series), GDS (GEO Datasets), GSM (GEO Samples) information for the NCBI GEO search. The GEO Platform record contains the list of elements on a gene expression microarray platform or elements that may be detected in a particular experiment. The GEO Samples record provides a detailed description for the experiment, such as experimental conditions and number of replicates, and must be associated with a GEO Platform record. The GEO Series record is a collection of all the samples in a given group enumerating the relationship and order (if any) of the samples to one another. The GEO Datasets record is a curated set of GEO Sample data and thus it may not be present for all experiments or organisms. GEO Datasets are associated based on the same platform and make it easy to directly compare platform elements, such as probe intensities. GEO_FTP_Addresses.txt is parsed by COnTORT into another file (GEO_FTP_Directories.txt) for use with the FTP retrieval process. COnTORT uses the information contained within these later two files to systematically retrieve all GDS, GSE, GSM, and GPL records associated with every search result provided. COnTORT exactly mimics the directory architecture of the NCBI FTP server. Note that COnTORT will not over-write GEO data already present at the destination.

After the retrieval of the NCBI GEO files is complete, COnTORT organizes the GEO Series data, which typically contains normalized transcriptomic data. This is accomplished by COnTORT running the GEOparse module to combine the annotations with the gene expression data file. All of the GEO Series files that were downloaded from NCBI GEO are copied to a new directory, unzipped, and then used to create the individual files that contain information for each experiment. Note that this file includes columns for every sample, even if the GEO Series record for that sample did not contain any gene expression data.

After the resulting GEOparse files are created, COnTORT correlates and combines the gene expression data based on gene identifiers present in the GenBank file. Using a custom Python script, COnTORT will identify all available gene annotations from the GenBank file (using the GenBank qualifiers “Locus Tag”, “Old Locus Tag”, “Gene Name”, and “Gene Synonyms”) in the GenBank file and use these to parse the transcriptomic data. Only those gene expression data sets associated with any of the annotation information present in the GenBank file are retained. This means that COnTORT will not include gene expression data that lacks clear annotation, including data organized by microarray probe ID or using an arbitrary gene designation ID that differs from those found in the GenBank file. While this can eliminate a small subset of gene expression data sets, it simplifies the downstream analysis of the remaining information. The transcript abundance data for each experiment is also mean normalized to correct for differences between different experiments and then all the files are concatenated into an COnTORT output that contains gene name/locus tag information. Note that blank entries are retained to indicate where gene expression data were missing in individual experiments.

The primary output of COnTORT are two single tab-delimited text files:  one file in which the data have been mean-centered for each experiment and one file in which the data have not been mean-normalized at all. Both files contain the results of the organization of the gene expression data. The first row contains the names of each experiment while the remaining rows contain the data for each gene present in the GenBank file. The first five columns contain the annotation information from the GenBank file (“Locus Tag”, “Old Locus Tag”, “Gene Name”, “Gene Synonyms”, and “Product”) with any values for a gene missing from the file identified as “N/A”. The remaining columns contain the gene expression data for each experimental sample. This output file can be opened in any text editor or spreadsheet program (such as Microsoft Excel) for manual processing or used in a variety of downstream applications programmatically.

Several additional files are created by COnTORT, and these are separated into separate new directories. All the NCBI GEO files retrieved in the search are located within the GEO Files directory. The GenBank annotations are organized in the ‘gene_annotations.txt’ file and located along wtih the FTP_download_log.txt file in the log_files directory. The GEOparse result files are located in the GEOannotate_results directory. Individual directories are also created to house files containing the matches between gene name/locus tag and expression data (GeneOrf_Match_Output) and the individual mean centered data files (Mean_Centered_Results). Finally files used to retrieve data from the NCBI FTP site are located with the FTP_Files directory. All files are retained locally for quality control and as controls for successful operation of COnTORT. Users can use these to trobuleshoot COnTORT if errors arrise.

### Tutorial for COnTORT

Running COnTORT is very user-friedly. There are two ways to obtain the program. The recommended way is to install it directly from PyPi using the command:

	pip install contort

This will allow you to call the program by opening the command line and simply typing:

	contort

Alternatively, you can clone the GitHub repository using the following command:

	git clone https://github.com/GLBRC/contort.git

Then to run COnTORT you can type the following into the command line (the location of contort.py will have to be indicated):

	python contort.py

Either way, COnTORT will present a Graphical User Interface (GUI) window for the uesr to provide the required information to run COnTORT. First the user will enter a search term that COnTORT will search at the NCBI GEO database and recover results for experiments involved in gene expression. This search will include all experiments categorized as high-throughput sequencing, expression profiling by array, and expression profiling by genome tiling array. Second, the user will open the GenBank file that will be used for annotating the results of the NCBI GEO search when combining the available gene expression data. The GenBank file can be downloaded from NCBI.

There are many GenBank files for every organism and so it is important to choose the appropriate GenBank annotation file for whatever search term they want to use. For example, there are many GenBank files for *Escherichia coli* for each genome that has been annotated. Generally users will use the standard GenBank annotation file for the reference genome, but other options are present at NCBI and it is important to choose the proper file so that the gene annotations from the NCBI GEO results will match the GenBank file.

After typing in the search term to use and selecting the GenBank file, COnTORT will run when the user selects the "Run" button. Additional information and progress is printed to the terminal window via standard out. At the end, COnTORT will report the length of time required to run.

To re-run COnTORT, you only need to re-run the program and enter a search term and select a GenBank file. It is recommendd that you create a new directory for every time you run COnTORT to avoid over-writing any files.
