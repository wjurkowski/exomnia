#!/usr/bin/perl -w
use strict;
use warnings;
use LWP;
use LWP::Simple;
use LWP::UserAgent;
my $Agent = LWP::UserAgent->new;

if ($#ARGV != 0) {die "Program used with parameters [input file]\n";}
#"Program used with parameters [fasta] [uniprot AC] [resn] [oAA] [nAA] [alignment] [tree]\n"

#paths
my $INST="/home/wiktor/Komoda/Projects/Aktualne/Epilepsy/netmol_pipeline/multiplex/";

#input data
my (@vpar,$fasta);
my @inputs=open_file($ARGV[0]);
foreach my $lin (@inputs){
 @vpar=split(/\t/,$lin);
 $fasta=$vpar[0];
 my $uniprotAC=$vpar[1];
 my $uniprotID=$vpar[2];
 my $resn=$vpar[3];
 my $oAA=$vpar[4];
 my $nAA=$vpar[5];
 my $alignment=$vpar[6];#in fasta format
 my $tree=$vpar[7];
 my $waryjat="$oAA$resn$nAA";# gv e.g. G124A

 #run Mutation Assesor
 #batch mode by web api
 #required input: uniport ID and variation
 #http://mutationassessor.org/v1/?cm=var&p=EGFR_HUMAN&var=G719S
 my (@MutAss);
 my $URL="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
 my $response = $Agent->post($URL);
 sleep 5;
 die "$URL error: ", $response->status_line
 unless  $response->is_success;
 my $output="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
 my $wyn = get($output);
 printf "$wyn\n"; 
 #print "$wyn\n";
 #@MutAss = split("\t",$wyn); 
 #printf "$MutAss[0]\t$MutAss[1]\n";

 #run PhD-SNP
 #python -O PhD-SNP.py -seq Test/1tthy.seq 21 K
 my $PhD="/usr/local/PhD-SNP2.0.6/";
 my $seq_simple=grep(!/>/, $fasta);
 `python -O $PhD/PhD-SNP.py -seq $seq_simple $resn $nAA >>PhD-SNP_output.txt`;

#run MAPP
#requires alignment in fasta format (ClustalW or ProbCons)
#requires tree build by ClustalW or Semphy
#output described in MAPP_readme.pdf
#physicochemical properties can be selected -s 1:3:4 (proterty 1+3+4), all on default 
#java -jar MAPP.jar -f LacI_Alignment.fa -t LacI.tree -o LacI_output.xls
#execute only if alignmnet and tree are present
 my $MAPP=$INST."MAPP";
 unless($alignment eq "-" and $tree eq "-"){
	`java -jar $MAPP/MAPP.jar -f $alignment -t $tree -o MAPP_output.txt`;
 }

#Create input files for others
 #Panther
 open (PTRIN, "> panther_csnpInput.txt") or die "Can not create Panther input file: $!"; 
 printf PTRIN "$waryjat\|$uniprotAC\|$resn\|$oAA\;$nAA\n";
 #PolyPhen
 open (PPHIN, "> pph_subs.input") or die "Can not create PPH input file: $!";
 printf PPHIN "$uniprotAC\t$resn\t$oAA\t$nAA\n"; 

}

#run PolyPhen
#PolyPhen-2 analysis pipeline consists of three separate components,
 #each one executed by a dedicated Perl program:
 #
 #  * MapSNPs     (mapsnps.pl)   Genomic SNP annotation tool
 #  * PolyPhen-2  (run_pph.pl)   Protein variant annotation tool
 #  * PolyPhen-2  (run_weka.pl)  Probabilistic variant classifier
 #$PPH/bin/run_pph.pl subs.pph.input 1>pph.features 2>run_pph.log &
 my $PPH=$INST."PolyPhen/polyphen-2.2.2";
 `$PPH/bin/run_pph.pl pph_subs.input 1>pph.features 2>run_pph.log &`;
 #$PPH/bin/run_weka.pl pph.features 1>pph.predictions
 `$PPH/bin/run_weka.pl pph.features 1>pph.predictions`;

#run Panther
#first classify protein sequence against PANTHER
 #pantherScore.pl -l <panther_hmm_library> -D B -V -i <fasta file> -o <output file> -n -T tmp/
 my $panther=$INST."Panther/csnpAnalysis1.02";
 `$panther/pantherScore.pl -l $panther/PANTHER7.2 -D B -V -i $fasta -o panther_scores.out -n -T tmp/`;
#classify SNPs
 #uprior.9comp contains some constants that could be modified
 #./snp_analysis.pl -l <panther_hmm_library> -c <score_outputfile_fromStepAbove> -s <csnp_input_file> -f <fasta_file> -b BLOSUM62 -V -p uprior.9comp  -o <output_file> -T tmp/ -a
 #If a SNP maps to multiple proteins, by default, the program will take the protein with the best HMm score.  Use the -a option if you want the results for all proteins
 #run analysis
print "$panther\n";
 `$panther/snp_analysis.pl -l $panther/PANTHER7.2 -c panther_scores.out -s panther_csnpInput.txt -f $fasta -b $panther/BLOSUM62 -V -p $panther/uprior.9comp  -o panther_snpanalysis.out -T tmp/`;
#output format in README
 my @outp=open_file("panther_snpanalysis.out");
 print "$outp[0]\t$outp[1]\t$outp[2]\t$outp[3]\t$outp[4]\t$outp[13]\t$outp[14]\t$outp[15]\t$outp[20]\n";

#other
sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}

