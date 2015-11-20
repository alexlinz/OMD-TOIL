###############################################################################
# readCountsCompete.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to count the MT reads which map to each gene in a
# set of reference genomes. Run this after you have run MTwrapperFunction.pl
# and readCounts.pl
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

################################################################################
### User-defined files and folder structure
################################################################################
my $genomeFolder = '../data/refGenomes';
my $gffFolder = '../data/gff';
my $mtFolder = '../data/AHC5MVBCXX_rRNA_processed';
my $mapFolder = '../data/mapping';
my $outputFolder = '../data/readCounts';

mkdir $outputFolder unless -d $outputFolder;

################################################################################
### Initialize lists of refGenomes and MTs
################################################################################

my @genomeList = glob($genomeFolder.'/*.fna');
for (@genomeList) {
  s/.+\/(.+).fna/$1/
  }

my @mtList;
opendir (DIR, $mtFolder);
while (my $dir = readdir(DIR)) {
  next if ($dir =~ m/^\./);
  push @mtList, $dir;
}
closedir(DIR);

################################################################################
### Generate files listing BAM files to aggregate
################################################################################

# Loop over each genome and MT
foreach my $mt (@mtList) {
# Create file  
  open(my $file, '>', $mt.'.txt');
  foreach my $genome (@genomeList) {
    # Write appropriate entry to file
    print $file $mapFolder."/".$mt."-".$genome.".sorted.bam\n";
  }
  close $file;
}

################################################################################
### Aggregate BAM files from the results of mapping
################################################################################

foreach my $mt (@mtList) {
  system("samtools merge -b ".$mt.".txt ".$mapFolder."/".$mt.".bam")
}

################################################################################
### Concatenate gff files
################################################################################

# Create output file
open(my $file, '>', $gffFolder.'/all.gff');

# Stream individual gff files to the all.gff file
foreach my $genome (@genomeList) {
  my $gff = $gffFolder."/".$genome.".gff";
  open(my $gffFile, '<', $gff);
  foreach my $line (<$gffFile>) {
    print $file $line;
  }
  close $gffFile;
}
close $file;

################################################################################
### Run htseq-count on the aggregated BAM files
################################################################################

foreach my $mt (@mtList) {
  system('/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t CDS -i locus_tag -m intersection-strict -o  '.$outputFolder.'/'.$mt.'.sam '.$mapFolder.'/'.$mt.'.bam '.$gffFolder.'/'.'all.gff >  '.$outputFolder.'/'.$mt.'.out');
}
