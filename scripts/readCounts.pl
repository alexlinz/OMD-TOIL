###############################################################################
# readCounts.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to count the MT reads which map to each gene in a
# genome. Run this after you have run MTwrapperFunction.pl
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
### Loop to obtain count data
################################################################################

# Loop over each genome and MT
foreach my $genome (@genomeList) {
    print "Processing genome ".$genome."\n";
    foreach my $mt (@mtList) {
      print "Processing metranscriptome ".$mt."\n";
      # Call HTseq-count
      system('/usr/local/bin/htseq-count -f bam -r pos -s no -a 0 -t CDS -i locus_tag -m intersection-strict -o '. $outputFolder.'/'.$mt.'-'.$genome.'.sam '.$mapFolder.'/'.$mt.'-'.$genome.'.sorted.bam '.$gffFolder.'/'.$genome.'.gff >  '.$outputFolder.'/'.$mt.'-'.$genome.'.out');
      }
  }
