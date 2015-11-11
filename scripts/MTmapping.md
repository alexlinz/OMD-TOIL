Mapping Metatranscriptomic Reads
===
Copyright (c) 2015, Joshua J. Hamilton  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Overview
--
Description of steps for mapping metatranscriptomic reads to our metagenomes.

Prerequisites
--
Ensure the following software is installed:  

* [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/) (BWA) - As of 2015-11-09. BWA version 0.7.5a-r405 is on Zissou and accessible by typing `bwa` in the command line.

* [samtools](http://www.htslib.org/download/) - This will need to be installed in your `home` folder.

Ensure the following data are available:  
* `pathToMTs` - a folder containing the metatranscriptomes to be mapped. As described in this repo's README, each sample should be in its own folder.

Create the following folders, if necessary  
* `pathToGenomes` -  a folder containing all reference genomes to be mapped. In FASTA nucleotide format.

Pipeline
--
This section describes the steps to map metatranscriptomic reads to the reference genomes. For convenience, the script `MTwrapperFunction.py` will execute the pipeline with a single command. The script takes as input the directories described below, and maps each metatranscriptome to each reference genome via the following commands:

1. Clean-up the FASTA header files. To avoid potential errors arising from non-alphanumeric characters (such as spaces or underscores), each contig is renamed to a number. For each genome `genome.fna` in `refGenomes`:

    `awk '/^>/{print ">" ++i; next}{print}' < genome.fna > genome.tmp`

    `mv genome.tmp genome.fna`

2. Index the reference genomes. For each genome `genome.fna` in `refGenomes`:

    `bwa index -p genome -a is genome.fna`

3. Map the metatranscriptomes to the reference genomes. For each metatranscriptome `sample_non_rRNA.fastq` and reference genome `genome.fna`:

    `bwa mem -t 30 genome.fna sample_non_rRNA.fastq > sample-genome.sam`

    where '-t' specifies the number of processors. Zissou has 32, please don't use all of them.

4. Manipulate the ouptut. Convert to BAM, sort, and index. Delete SAM and unsorted BAM files to save space. For each `sample-genome.sam` file:

    `samtools view -b  -S -o sample-genome.bam sample-genome.sam`

    `samtools sort -o sample-genome.sorted.bam -O bam -T /temp/sample-genome sample-genome.bam`

    `samtools index sample-genome.sorted.bam`

    `rm sample.sam`

    `rm sample.bam`

The script is called as follows:

    `perl MTwrapperFunction pathToMTs pathToGenomes pathToMaps numProcs`

The inputs to the wrapper function are as follows:

  | Argument | Description  |
  |---|---|
  | pathToMTs | location of the metatranscriptome reads to be mapped |
  | pathToGenomes | location of the reference genomes |
  | pathToMaps | desired output location of the mapped reads |
  | numProcs | number of processors to use for mapping |
