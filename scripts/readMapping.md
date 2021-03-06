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

* [Genome Tools](http://genometools.org/pub/) - This will need to be installed in your `home` folder.

* [samtools](http://www.htslib.org/download/) - This will need to be installed in your `home` folder.

* [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html) - As of 2015-11-17, this Python package is installed on Zissou.  We will be using the script `htseq-count` which is accessible from the command line.

Ensure the following data are available:  

* `pathToMTs` - a folder containing the metatranscriptomes to be mapped. As described in this repo's README, each sample should be in its own folder.

Create the following folders, if necessary  

* `pathToGenomes` -  a folder containing all reference genomes to be mapped. In FASTA nucleotide format (`file.fna`). __Note__: The FASTA nucleotide sequence for the internal standard should also be included here.

* `pathToGFFs` -  a folder containing all GFF files of all reference genomes to be mapped. In FASTA nucleotide format. __Note__: A dummy GFF file for the internal standard should also be included here.

Data Processing
---
In order for `htseq-count` to properly count reads, the GFF files must be in the proper format. The script `gffValidator` will analyze all GFF files in the `pathToGFFs` folder, and generate a file `pathToGFFs/all.out` describing errors in the GFF files. In Josh's experience (thus far), errors will be one of the following messages:

* `token "string" on line XYZ in file “genome.gff" does not contain exactly one ‘=‘`

  This error occurs when a gene annotation contains a semi-colon. According to the GFFv3 [standard](http://www.sequenceontology.org/gff3.shtml), semi-colons have a reserved meaning and must be escaped. The standard recommends using the URL escape `%3B`.

* `error: could not parse score '' on line XYZ in file ‘genome.gff’`

  This error occurs because a line of the GFF file is incomplete (e.g., missing the score information). Josh found this to occur only for CRISPR arrays defined in the final line of the GFF file. Since he isn't concerned with mapping to CRISPR arrays, he deleted these lines.

Mapping of Reads
--
This section describes the steps to map metatranscriptomic reads to the reference genomes. For convenience, the script `MTwrapperFunction.pl` will execute the pipeline with a single command. The script takes as input the directories described below, and maps each metatranscriptome to each reference genome via the following commands:

1. Index the reference genomes. For each genome `genome.fna` in `refGenomes`:

    `bwa index -p genome -a is genome.fna`

2. Map the metatranscriptomes to the reference genomes. For each metatranscriptome `sample_non_rRNA.fastq` and reference genome `genome.fna`:

    `bwa mem -t 30 genome.fna sample_non_rRNA.fastq > sample-genome.sam`

    where '-t' specifies the number of processors. Zissou has 32, please don't use all of them.

3. Manipulate the ouptut. Convert to BAM, sort, and index. Delete SAM and unsorted BAM files to save space. For each `sample-genome.sam` file:

    `samtools view -b  -S -o sample-genome.bam sample-genome.sam`

    `samtools sort -o sample-genome.sorted.bam -O bam -T /temp/sample-genome sample-genome.bam`

    `samtools index sample-genome.sorted.bam`

    `rm sample.sam`

    `rm sample.bam`

The script is called as follows:

    `perl MTwrapperFunction pathToMTs pathToGenomes pathToMaps numProcs`

__Note:__ For convenience, the script is currently hard-coded to use the folder structure described in this repo. The script also specifies use of 24 processors.

The inputs to the wrapper function are as follows:

| Argument | Description  |
|---|---|
| pathToMTs | location of the metatranscriptome reads to be mapped |
| pathToGenomes | location of the reference genomes |
| pathToMaps | desired output location of the mapped reads |
| numProcs | number of processors to use for mapping |

Counting of Mapped Reads: Uncompetitive
--

This section describes the steps to count the metatranscriptomic reads which mapped to each gene in a reference genome. __Because each genome is processed individually, it is possible that some reads are being counted in multiple genomes.__ For convenience, the script `readCounts.pl` will execute the pipeline with a single command. The script takes as input the directories described below, and counts reads mapped to each gene using [htseq-count](http://www-huber.embl.de/HTSeq/doc/count.html#count) for each (metatranscriptome, genome) pair:

    `/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t CDS -i locus_tag -m intersection-strict -o outputFolder/mt-genome.sam mapFolder/mt-genome.sorted.bam gffFolder/genome.gff > outputFolder/mt-genome.out`

where variables are defined as follows:

| Argument | Description  |
|---|---|
| outputFolder | location to store read counts for each genome |
| gffFolder | location of the gff files for the reference genomes |
| mapFolder | location of the mapped reads |
| mt | the current metatranscriptome|
| genome | the current genome |

and flags are defined as follows:

| Flag | Description  |
|---|---|
| -f bam | mapped reads are stored in a `bam` file |
| -r pos | `bam` file has been sorted by position along the genome |
| -s no | reads are not strand-specific (e.g., they are from cDNA and can match to the same or opposite strand as the gene) |
| -t CDS | type of feature to look for, options available in the the `gff` file |
| -i locus_tag | label to apply to each gene, options available in the the `gff` file |
| -m intersection-strict | mapped reads must fully lie within a single gene |
| -o | location of output `sam` file|
| > | location of file with read counts |

Additional info about these flags is available in the `htseq-count` [documentation](http://www-huber.embl.de/HTSeq/doc/count.html#count) and the `GFF` [file specification](http://gmod.org/wiki/GFF2).

The output of the script is a tab-delimited `.out` file giving the locus tag and number of mapped reads.

The script is called as follows:

    `perl readCounts.pl genomeFolder gffFolder mtFolder mapFolder outputFolder`

The inputs to the wrapper function are as follows:

| Argument | Description  |
|---|---|
| genomeFolder | location of the reference genomes |
| gffFolder | location of the gff files for the reference genomes |
| mtFolder | location of the metatranscriptome reads |
| mapFolder | location of the mapped reads |
| outputFolder | location to store read counts for each genome |

__Note:__ For convenience, the script is currently hard-coded to use the folder structure described in this repo.

Counting of Mapped Reads: Competitive
--

In the previous section, `htseq-count` was run on each individual (metatranscriptome, genome) pair. As a result, it is possible that some reads mapped to genes in different genomes. To perform a competitive mapping (e.g., count only those reads which map unambiguously to a single gene), the `bam` files for all (metatranscriptome, genome) pairs from the same metatranscriptome can be aggregated, and `htseq-count` run on each aggregated set.

_This analysis needs to be performed and written up._

Counting of Mapped Reads
--

  The python script `processReadCounts.py` aggregates these results into a single `.csv` file, showing the percent of reads in each MT which mapped to each genome.

_This analysis needs to be updated with proper normalization based on the internal standard._
