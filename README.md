OMD-TOIL
===
Copyright (c) 2015, Katherine D. McMahon and Lackeys  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Overview
--
During August 20-21, 2015, the McMahon Lab carried out a 24 hour sampling of Lake Mendota. This repo describes the data and analysis associated with this expedition, known as OMD-TOIL (Operation Mendota Drain: Transcriptomics of Inland Lakes).

Samples were collected approximately every two hours from 6am on August 20 to 6am on August 21. A subset of samples were subject to RNA extraction, spiking with an internal standard, and sent to the [UW Biotech Center](https://www.biotech.wisc.edu/) for sequencing. Resulting sequences were subject to rRNA removal.

Repo Structure
--
    project
    |- README                                             # the top level description of content
    |- data                                               # organizational structure of raw and primary data
    | |- AHC5MVBCXX                                       # raw sequence reads
    | | |- McMahon_Demultiplex_Stats.htm                  # QC information for each sample
    | | |- ME150256                                       # raw sequence reads for sample ME150256
    | | | |- ME15025616A_ACAGTG_L001_R1_001_fastqc.html   
    | | | |- ME15025616A_ACAGTG_L001_R1_001.fastq.gz      
    | | | |- ME15025616A_ACAGTG_L001_R1_001.fastq.gz.md5  
    | | | |- ME15025616A_ACAGTG_L001_R2_001_fastqc.html   
    | | | |- ME15025616A_ACAGTG_L001_R2_001.fastq.gz      
    | | | |- ME15025616A_ACAGTG_L001_R2_001.fastq.gz.md5  
    | | |- ME150270                                       # raw sequence reads for sample ME150270. Same file structure as for ME150256.
    | | |- ME150276                                       # raw sequence reads for sample ME150276. Same file structure as for ME150256.
    | | |- ME150286                                       # raw sequence reads for sample ME150286. Same file structure as for ME150256.
    | |- AHC5MVBCXX_merged                                # merged sequence reads
    | | |- ME150256                                       # merged sequence reads for sample ME150256
    | | | |- out.extendedFrags.fastq
    | | | |- out.hist
    | | | |- out.histogram
    | | | |- out.notCombined.fastq
    | | |- ME150270                                       # merged sequence reads for sample ME150270. Same file structure as for ME150256.
    | | |- ME150276                                       # merged sequence reads for sample ME150276. Same file structure as for ME150256.
    | | |- ME150286                                       # merged sequence reads for sample ME150286. Same file structure as for ME150256.
    | |- AHC5MVBCXX_rRNA_processed                        # sequence reads with rRNA removed
    | | |- ME150256                                       # sequence reads with rRNA removed for sample ME150256
    | | | |- ME15025616A_non_rRNA.fastq  
    | | | |- ME15025616A_rRNA.log  
    | | | |- out.extendedFrags.fastq  
    | | | |- out.histogram
    | | | |- ME15025616A_rRNA.fastq      
    | | | |- ME15025616A_rRNA.sam  
    | | | |- out.hist                 
    | | | |- out.notCombined.fastq
    | | |- ME150270                                       # sequence reads with rRNA removed for sample ME150270. Same file structure as for ME150256.
    | | |- ME150276                                       # sequence reads with rRNA removed for sample ME150276. Same file structure as for ME150256.
    | | |- ME150286                                       # sequence reads with rRNA removed for sample ME150286. Same file structure as for ME150256.
    |- metadata                                           # metadata associated with the samples
    |- protocols                                          # protocols associated with sample collection and processing
    | |- SamplingProcedures.pdf                           # description of sampling procedures
    | |- RNAExtraction.txt                                # RNA extraction for sequencing
    | |- Sequencing.txt                                   # Library prep and sequencing protocols, if available. Otherwise, a description of methods.
    |- scripts                                            # protocols associated with processing and analysis of sequence data

Notes on Data Storage
--
The sequencing data associated with this project are quite large and not stored in this repo. As of 2015-11-09, these dare are stored on Zissou at `/home/fmoya/RNAseq_OMD_TOIL15`.
