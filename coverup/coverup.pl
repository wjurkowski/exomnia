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
	  if($sample{$var} >= $xcov){
		$knownHQ{$zn} = $sample{$zn};
		print #print known HQ
	  }
	  else{
		$knownLQ{$zn} = $sample{$zn};
		print #print known LQ
	  }
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
	  @=split("-",$var);
      my $pos=~ /\d+/;
      my $eb=$start{$gene}{$exon};
      my $ee=$stop{$gene}{$exon};
      if($pos >= $eb and $pos <= $ee){
		if($qual[1] = 100.00){#variant could be potentially OK - store confirmed missing known epi variants
			$nima{$var}=$equal{$gene}{$exon};
		}
		else{#low quality data
			$couldbe{$var}=$equal{$gene}{$exon};
			#covered but low quality
			$readd=get_read_depth();
		}
      }
    }
  }
  printf OUT1 "$gene\t$gq\n";#print short summary of genes quality
}

sub get_read_depth {
my $key=shift;
`samtools faidx /work/projects/cogie/analysis/func2vcf-tests/reference/hg19.fa chr1:2000001-2000002`;
}

sub ecover {#read in exon coverages
my @exoncover=shift;
my (%start, %stop, %minDepth, %avDepth, %maxDepth, %exonCover);
foreach $variant(@exoncover){
  my @data=split("\t",$variant);
  $symbol=$data[0];
  $exon=$data[1];
  $chr=$data[2];
  $enst=$data[5];
  $ensg=$data[6];
  $gene=$symbol."-".$chr."-".$ensg."-".$enst;
  $start{$gene}{$exon}=$data[3];
  $stop{$gene}{$exon}=$data[4];
  $minDepth=$data[7];
  $avgDepth=$data[8];
  $maxDepth=$data[9];
  $exonCover=$data[10];
  
  #for each exon calculate read depth-based score:
  #	s=(max-min)/6
  #	score=|avg-((max-min)/6)|
  #	threshold of the score estimated based on normal distribution of read depth. If > ~95% population have read depth >= 10 then average should be around 30.
  #	score for such a distribution (min=0, max=60, avg=30) is 10
  $s=($maxDepth-$minDepth)/6;
  $score=$avgDepth-2*$s;
  if($score > 10 and $exonCover > 95.00){$qcat=1.0;}
  if($score > 10 and $exonCover > 90.00){$qcat=0.8;}
  if($score > 10 and $exonCover > 85.00){$qcat=0.6;}
  if($score > 10 and $exonCover > 50.00){$qcat=0.4;}
  if($score > 5 and $exonCover > 30.00){$qcat=0.2;}
  if($exonCover < 50.00){$qcat=0.0;}
  $qa=$score.",".$exonCover.",".$qcat;
  $equal{$gene}{$exon}=$qa;
#exons quality is defined in following classes: 
#	very high: all exons with of high quality ;
#	high: all exons with at least high quality
}
return \%equal,\%start,%stop;
#end of equal function
}

sub readfunc{#read in func file
my @dataset=shift;
my (%sample);
foreach my $variant(@dataset){#for each data set
    unless($variant =~ /^#/){
      my @data=split("\t",$variant);
      my $chrom=$data[0];
      my $pos=$data[1];
      my $ref=$data[3];
      my $alt=$data[11];
      #correct empty field
      if($ref eq "-"){$ref="null";}
      if($alt eq "-"){$alt="null";}
      #correct homopolymer indels
      if($ref=~/\(/){
		$ref=~s/\(//;
		$ref=~s/\)/;/;
		my @hind=split(";",$ref);
		$ref=~s/$hind[0]\;\d+/$hind[0] x $hind[1]/e;
      }
      if($alt=~/\(/){
		$alt=~s/\(//;
		$alt=~s/\)/;/;
		my @hind=split(";",$alt);
		$alt=~s/$hind[0]\;\d+/$hind[0] x $hind[1]/e;
      }
      my $type=$data[5];
      my $ensg=$data[21];
      my $enst=$data[22];
      my $ensp=$data[23];
      my $hgnc=$data[26];
      my $cdna=$data[34];
      my $prot=$data[35];
	  my $var=$hgnc."-".$ref.$pos.$alt; #each variation is a key
      #
      #skip ROH variants
      if($type eq "ROH"){next;}
      #correct incorrect delimiter in genotype field
      $data[13]=~s/:/\//;#replace wrong allele delimiter
      $data[13]="./." if $data[13] eq "";
      #get allele counts	
      my $sAN=2;
      my $sAC=0;
      if($data[13] eq "1|1" or $data[13] eq "1/1"){$sAC=2;}
      elsif($data[13] eq "0|1" or $data[13] eq "0/1"){$sAC=1;}
	  #non-redundant variants only 	
      unless($genotype{$var}){
		#derivatives
		$AN{$var}=0 unless $AN{$var};
		$AC{$var}=0 unless $AC{$var};
		$NS{$var}++;#total number of samples/data sets
		$AN{$var}=$AN{$var}+$sAN;#total number of alleles
		$AC{$var}=$AC{$var}+$sAC;#total number allele counts
		#genotype
		$genotype{$var}=$data[13];
		my $cover=$data[8];
		my $gtpfield=$genotype{$var}.":".$cover;
		#define hash of arrays to store genotype info
		#add new array value for each func file
		push (@{$genotyping{$var}},$gtpfield);
	  }
	  unless($sample{$var}){$sample{$var}=$data[8];}
    }
}
return %sample;
}

#read in known epilepsy-associated variants
sub epigv{
my @epi=shift;
foreach my $variant(@epi){#for each data set
	my @data=split("\t",$variant);
	$known{$data[0]}=$data[3];#uniprot AC
}
return %known;
}

#read files
sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}

