###############################################################################
# MTwrapperFunction.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to map metagenomics reads to reference genomes
# using BWA. It was developed to map reads from OMD-TOIL to our SAGs and MAGs.
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

################################################################################
### User-defined files and folder structure
################################################################################
#my $usage  = "Command sequence: perl MTwrapperFunction.pl pathToMTs pathToGenomes pathToMaps numProcs \n";
#my $pathToMTs  = shift or die $usage;
#my $pathToGenomes  = shift or die $usage;
#my $pathToMaps  = shift or die $usage;
#my $numProcs  = shift or die $usage;

my $pathToMTs = '../data/AHC5MVBCXX_rRNA_processed';
my $pathToGenomes = '../data/refGenomes';
my $pathToMaps = '../data/mapping';
my $numProcs = 24;

my @genomeList = glob($pathToGenomes.'/*.fna');
my @sampleList = glob($pathToMTs.'/*');

################################################################################
### Step 1: Clean up FASTA Header Sequences
################################################################################

print "Renaming FASTA header lines.\n";

foreach my $genomePath (@genomeList) {
  system("awk '/^>/{print \">\" ++i; next}{print}' < ".$genomePath." >".$genomePath.".proc");
  system("mv ".$genomePath.".proc ".$genomePath);
}

################################################################################
### Step 2: Index Genomes
################################################################################

my $int = 1;
foreach my $genomePath (@genomeList) {
  print "Indexing genome ".$int." of ".@genomeList."\n";
  if ($genomePath =~ /(.+).fna/) {
    my $genome = $1;
    system("bwa index -p $genome -a is $genomePath");
    $int++;
    }
}

################################################################################
### Step 3: Map Metatranscriptomes to Reference Genomes
################################################################################

$int = 1;
my $total = @genomeList*@sampleList;

foreach my $samplePath (@sampleList) {
  foreach my $genomePath (@genomeList) {
    if ($genomePath =~ /.+\/(.+).fna/) {
      my $genome = $1;
      if ($samplePath =~ /.+\/(.+)/) {
	my $sample = $1;
	print "Mapping sample ".$sample." against genome ".$genome." (".$int." of ".$total."). \n";
	$int ++;
	system("bwa mem -t ".$numProcs." ".$pathToGenomes."/".$genome." ".$samplePath."/".$sample."_non_rRNA.fastq > ".$pathToMaps."/".$sample."-".$genome.".sam");
      }
    }
  }
}

################################################################################
### Step 4: Manipulate Output SAM files. Convert to BAM, sort and index. Delete
### unsorted SAM and BAM files to save space.
################################################################################

my @samList = glob($pathToMaps.'/*.sam');
$int = 1;

foreach my $samPath (@samList) {
  print "Processing SAM file ".$int." of ".@samList."\n";
  if ($samPath =~ /.+\/(.+).sam/) {
    my $sam = $1;
    system("samtools view -b -S -o ".$pathToMaps."/".$sam.".bam ".$pathToMaps."/".$sam.".sam");
    system("samtools sort -o ".$pathToMaps."/".$sam.".sorted.bam -O bam -T ".$pathToMaps."/temp/".$sam." ".$pathToMaps."/".$sam.".bam");
    system("samtools index ".$pathToMaps."/".$sam.".sorted.bam");
    system("rm ".$pathToMaps."/".$sam.".bam");
    system("rm ".$pathToMaps."/".$sam.".sam");
    $int ++;
  }
}


