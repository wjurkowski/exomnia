#!/usr/bin/perl -w
use strict;
use warnings;
use LWP;
use LWP::Simple;
use LWP::UserAgent;
my $Agent = LWP::UserAgent->new;

if ($#ARGV != 0) {die "Program used with parameters [input file]\n";}

#paths
my $INST="/home/wiktor/Komoda/Projects/Aktualne/Epilepsy/netmol_pipeline/multiplex/";
#
#input data
#my (@vpar,$fasta);
my @inputs=open_file($ARGV[0]);
foreach my $lin (@inputs){
 @vpar=split(/\t/,$lin);
 $gb=$vpar[0];#genome build, hg19 default
 my $chr=$vpar[1];#chromosome
 my $pos=$vpar[2];#position
 my $basen=$vpar[3];
 my $oBP=$vpar[4];
 my $nBP=$vpar[5];
 my $waryjat="$chr,$pos,$oBP,$nBP";# gv e.g. G124A

#run Mutation Assesor
#batch mode by web api
#required input: chromosome and position and variation
#http://mutationassessor.org/v1/?cm=var&var=7,55178574,G,A
#http://mutationassessor.org/v1/?cm=var&var=hg19,13,32912555,G,T 
 my (@MutAss);
 my $URL=`wget http://mutationassessor.org/v1/?cm=var&var=$waryjat&frm=txt`;
 my $response = $Agent->post($URL);
 sleep 5;
 die "$URL error: ", $response->status_line
 unless  $response->is_success;
 my $wyn = get($output);
 printf "$wyn\n";


#run MutationTaster
#batch mode by web api
#required input: file with snippets
#perl mutation_taster_batch_query.pl -i snippets_example.tsv -o outdir
my $snippets_input=;
#It can be generated from VCF files using other available perl script 
#ShrinkMasterVarFile.pl -i masterVarBeta-HG00732-200-37-ASM.tsv -homozygous_only true -mincov 40
`perl mutation_taster_batch_query.pl -i $snippets_input -o MutationTaster_outdir`;
#analyze results
#perl mutation_taster_results.pl [-d dir]
`perl mutation_taster_results.pl -d outdir`;
#results

#run snpEff
#requires input file in the VCF format or in the text format
#java -jar snpEff.jar -c snpEff.config hg19 -i vcf  example_vcf.txt
#java -jar snpEff.jar -c snpEff.config hg19 -i txt  example_vcf.txt
#txt format: 5   140532    T    C 
`java -jar snpEff.jar -c snpEff.config hg19 snps.txt`;

