#!/usr/bin/perl -w
use strict;
use warnings;

if ($#ARGV != 3) {die "Program used with parameters [fasta] [uniprot AC] [resn] [oAA] [nAA] [alignment] [tree]\n";}

#input data
my $fasta=$ARGV[0];
my $uniprotAC=$ARGV[1];
my $resn=$ARGV[2];
my $oAA=$ARGV[3];
my $nAA=$ARGV[4];
my $alignment=$ARGV[5];#in fasta format
my $tree=$ARGV[6];
my $waryjat="$oAA.$resn.$nAA";# gv e.g. G124A

#run Mutation Assesor
#batch mode by web api
#required input: uniport ID and variation
#http://mutationassessor.org/v1/?cm=var&p=EGFR_HUMAN&var=G719S
my (@MutAss);
my $wyn=`wget http://mutationassessor.org/v1/?cm=var&p=$uniprotAC&var=$waryjat&frm=txt`;
print "$wyn\n";
@MutAss = split("\t",$wyn); 
printf "$MutAss[0]\t$MutAss[1]\n";

#run PolyPhen
#PolyPhen-2 analysis pipeline consists of three separate components,
#each one executed by a dedicated Perl program:
#
#  * MapSNPs     (mapsnps.pl)   Genomic SNP annotation tool
#  * PolyPhen-2  (run_pph.pl)   Protein variant annotation tool
#  * PolyPhen-2  (run_weka.pl)  Probabilistic variant classifier
#$PPH/bin/run_pph.pl subs.pph.input 1>pph.features 2>run_pph.log &
my $PPH="/home/wiktor/Komoda/Projects/Aktualne/Epilepsy/netmol_pipeline/multiplex/PolyPhen/polyphen-2.2.2";
open (PPHIN, "> pph_subs.input") or die "Can not create PPH input file: $!";
printf PPHIN "$uniprotAC\t$resn\t$oAA\t$nAA\n"; 
`$PPH/bin/run_pph.pl pph_subs.input 1>pph.features 2>run_pph.log &`;
#$PPH/bin/run_weka.pl pph.features 1>pph.predictions
`$PPH/bin/run_weka.pl pph.features 1>pph.predictions`;

#run Panther
#first classify protein sequence against PANTHER
 #pantherScore.pl -l <panther_hmm_library> -D B -V -i <fasta file> -o <output file> -n -T tmp/
`pantherScore.pl -l PANTHER7.2 -D B -V -i $fasta -o panther_scores.out -n -T tmp/`;
#classify SNPs
 #uprior.9comp contains some constants that could be modified
 #./snp_analysis.pl -l <panther_hmm_library> -c <score_outputfile_fromStepAbove> -s <csnp_input_file> -f <fasta_file> -b BLOSUM62 -V -p uprior.9comp  -o <output_file> -T tmp/ -a
 #If a SNP maps to multiple proteins, by default, the program will take the protein with the best HMm score.  Use the -a option if you want the results for all proteins
 #create input file
 open (PTRIN, "> panther_csnpInput.txt") or die "Can not create Panther input file: $!"; 
 printf PTRIN "$waryjat\|$uniprotAC\|$resn\|$oAA\;$nAA\n";
 #run analysis
 `snp_analysis.pl -l PANTHER7.2 -c panther_scores.out -s panther_csnpInput.txt -f test.fasta -b BLOSUM62 -V -p uprior.9comp  -o panther_snpanalysis.out -T tmp/`;
#output format in README

#run MAPP
#requires alignment in fasta format (ClustalW or ProbCons)
#requires tree build by ClustalW or Semphy
#output described in MAPP_readme.pdf
#physicochemical properties can be selected -s 1:3:4 (proterty 1+3+4), all on default 
#java -jar MAPP.jar -f LacI_Alignment.fa -t LacI.tree -o LacI_output.xls
#execute only if alignmnet and tree are present
unless($alignment eq "-" and $tree eq "-"){
	`java -jar MAPP.jar -f $alignment -t $tree -o MAPP_output.txt`;
}

#run PhD-SNP
#python -O PhD-SNP.py -seq Test/1tthy.seq 21 K
my $seq_simple=grep(!/>/, $fasta);
`python -O PhD-SNP.py -seq $seq_simple $resn $nAA`;


