# Script to parse vcf file and save variants with allele frequency higher than specified
# 
# Copyright:: Copyright 2012-2013 University of Luxembourg/Luxembourg Centre for Systems Biomedicine
# Authors::   Wiktor Jurkowski <wiktor.jurkowski@uni.lu>

#!/usr/bin/perl -w
use strict;
use warnings;
use Math::Round qw(nearest);

if ($#ARGV != 2) {die "Program requires parameters! [exon coevrage] [func file] [known epi] [coverage threshold]\n";}
#per individual: what is the fraction of epi snps that can't be potentially identified because low exon coverage? 
#for each data set print for each gene: table with columns: symbol	exon_quality

#read in exon coverage file
my @exons=open_file($ARGV[0]);
#read in corresponding func file
my @func=open_file($ARGV[1]);
#read in known epilepsy variants in exomnia format
my @epi=open_file($ARGV[2]);

#read in options
my $xcov=$ARGV[3];

#define output files
open (OUT1, ">highQ_genes.txt"); 
open (OUT2, ">.txt"); 

my %sample = readfunc(@func);#hash of variants, key -> key, value ->  coverage (read depth)
my %known = epigv(@epi);#hash of variants, key -> key, value ->  gv type
my (%equal,%start,%stop) = ecover(@exons);#hashes (of hashes) for each gene and each exon storing quality info, start and stop chr positions

#check coverage of known variants  
#
#check if known variant is in hash
for my $var(keys %sample){
  for my $zn(keys %known){
    if(%sample{$zn}){#if present check read depth
	  if($sample{$var} >= $xcov){$knownHQ{$zn} = $sample{$zn};}
	  else{$knownLQ{$zn} = $sample{$zn};}
    }
    else{$knownmissing{$zn}=$known{$zn};}
  }
}
#
#check coverage of that particular exon to see if any variant from this region could be potentially included
#
for my $gene (keys %equal){
  my $gq="high";
  for my $exon (keys %{ $equal{$gene}} ){
    my @qual = split(",",$equal{$gene}{$exon});
    if($qual[2] < 0.6){$gq="low";}#check if all exons have reasonable quality
    #check known missing variants
    for my $var {keys %knownmissing}{
	  if($readdepth >= $xcover){#any other variant at that position?
		$realref=get_reference();	
	  }
      if($pos >= $eb and $pos <= $ee){
		if($qual[1] >= 85.00){#variant could be potentially OK - store confirmed missing known epi variants
			$nima{$var}=$equal{$gene}{$exon};
		}
		else{#low quality data
			$couldbe{$var}=$equal{$gene}{$exon};
			#covered but low quality
			#$readd=get_read_depth();
		}
      }
    }
  }
  printf OUT1 "$gene\t$gq\n";#print short summary of genes quality
}


