OMD-TOIL
===
Copyright (c) 2015, Katherine D. McMahon and Lackeys  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Overview
--
During August 20-21, 2015, the McMahon Lab carried out a 24 hour sampling of Lake Mendota. This repo describes the data and analysis associated with this expedition, known as OMD-TOIL (Operation Mendota Drain: Transcriptomics of Inland Lakes).

Samples were collected approximately every two hours from 6am on August 20 to 6am on August 21. A subset of samples were subject to RNA extraction, spiking with an internal standard, and sent to the [UW Biotech Center](https://www.biotech.wisc.edu/) for sequencing. Resulting sequences were subject to rRNA removal. These sequences were then mapped against our collection of reference genomes using BWA.

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
    | |- AHC5MVBCXX_merged                                # merged sequence reads
    | | |- ME150256                                       # merged sequence reads for sample ME150256
    | | | |- out.extendedFrags.fastq
    | | | |- out.hist
    | | | |- out.histogram
    | | | |- out.notCombined.fastq
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
    | |- mapping                                          # results of mapping metatranscriptomes against reference genomes
    | | |- ME15025616A-2236347068.sorted.bam              # sorted BAM file of MT sample ME15025616A mapped against genome 2236347068
    | | |- ME15025616A-2236347068.sorted.bam.bai          # index for the above BAM file
    | |- refGenomes                                       # reference genomes against which metatranscriptomes will get mapped
    | | |- 2236347068.fna                                 # unannotated genome 2236347068. The number 2236347068 is the IMG OID.
    | | |- 2236347068.amb                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.ann                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.bwt                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.pac                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.sa                                  # file generated by bwa in preparation for mapping
    | | |- 2236347068.genes.fna                                 # annotated genome 2236347068. The number 2236347068 is the IMG OID.
    | | |- 2236347068.genes.amb                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.genes.ann                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.genes.bwt                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.genes.pac                                 # file generated by bwa in preparation for mapping
    | | |- 2236347068.genes.sa                                  # file generated by bwa in preparation for mapping
    |- metadata                                           # metadata associated with the samples
    |- protocols                                          # protocols associated with sample collection and processing
    | |- SamplingProcedures.pdf                           # description of sampling procedures
    | |- RNAExtraction.txt                                # RNA extraction for sequencing
    | |- Sequencing.txt                                   # library prep and sequencing protocols, if available. Otherwise, a description of methods.
    |- scripts                                            # protocols associated with processing and analysis of sequence data
    | |- merging.txt                                      # script/command for merging raw reads.
    | |- MTmapping.md                                     # workflow for mapping metatranscriptomic reads
    | |- MTwrapperFunction.py                             # python wrapper function which executes most steps in the MTmapping.md workflow
    | |- remove_rRNA.txt                                  # script/command for removing rRNA sequences.
Notes on Data Storage
--
The sequencing data associated with this project are quite large and not stored in this repo. As of 2015-11-09, these dare are stored on Zissou at `/home/fmoya/RNAseq_OMD_TOIL15`.
