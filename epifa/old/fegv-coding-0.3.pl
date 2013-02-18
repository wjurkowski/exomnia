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

#output file
open (GOUT, ">SNPs_function_pipeline-out.txt") or die "Can not create general output file: $!";

#-----------------------create batch input for PPH and Panther---------------
#-----------------------run remaining tools-----------------------------------
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
 open (MAOUT, "> MutationAssesor_output.txt") or die "Can not create MutAss output file: $!";
 my $URL="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
 my $response = $Agent->post($URL);
 sleep 5;
 die "$URL error: ", $response->status_line
 unless  $response->is_success;
 my $output="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
 my $wyn = get($output);
 printf MAOUT "$wyn\n"; 

 #run PhD-SNP
 #python -O PhD-SNP.py -seq Test/1tthy.seq 21 K
 open (PDOUT, "> PhD-SNP_output.txt") or die "Can not create PhD-SNP output file: $!";
 my $PhD="/usr/local/PhD-SNP2.0.6/";
 my $seq_simple=grep(!/>/, $fasta);
 `python -O $PhD/PhD-SNP.py -seq $seq_simple $resn $nAA > outtmp`;
 my @out=open_file("outtmp");
 my @PhDout=split(/ \+/,$out[10]);
 print "$PhDout[4]\n";
 printf PDOUT "$PhDout[0]\t$PhDout[1]\t$PhDout[2]\t$PhDout[3]\t$PhDout[4]\n"; 

#run MAPP
#requires alignment in fasta format (ClustalW or ProbCons)
#requires tree build by ClustalW or Semphy
#output described in MAPP_readme.pdf
#physicochemical properties can be selected -s 1:3:4 (proterty 1+3+4), all on default 
#java -jar MAPP.jar -f LacI_Alignment.fa -t LacI.tree -o LacI_output.xls
#execute only if alignmnet and tree are present
 my $MAPP=$INST."MAPP";
 unless($alignment eq "-" and $tree eq "-"){
   `java -jar $MAPP/MAPP.jar -f $alignment -t $tree -o MAPP_tmp.txt`;
   `cat MAPP_tmp.txt >> MAPP_output.txt`;
 }

#Create input files for others
 #Panther
 open (PTRIN, "> panther_csnpInput.txt") or die "Can not create Panther input file: $!"; 
 printf PTRIN "$waryjat\|$uniprotAC\|$resn\|$oAA\;$nAA\n";
 #PolyPhen
 open (PPHIN, "> pph_subs.input") or die "Can not create PPH input file: $!";
 printf PPHIN "$uniprotAC\t$resn\t$oAA\t$nAA\n"; 

}

#------------------------------run Panther and PPH------------------------
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
 my @PPHout=open_file("pph.predictions");
 #output legend
 #http://genetics.bwh.harvard.edu/pph2/dokuwiki/appendix_a
 #dbSNP prediction pph2_class prob FPR sensitivity FDR dScore Score1 Score2 ident lenght normSAS SS MapReg dVol dProp Transv CpG MinDJnc Pfam
 printf GOUT "$PPHout[4]\t$PPHout[11]\t$PPHout[14]\t$PPHout[15]\t$PPHout[16]\t$PPHout[17]\t$PPHout[18]\t$PPHout[22]\t$PPHout[23]\t$PPHout[24]\t$PPHout[32]\t$PPHout[33]\t$PPHout[34]\t$PPHout[35]\t$PPHout[36]\t$PPHout[37]\t$PPHout[38]\t$PPHout[47]\t$PPHout[49]\t$PPHout[50]\t$PPHout[51]\t";


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
 my @PTRout=open_file("panther_snpanalysis.out");
 #output legend
 #snpId seqId subPSEC Pdeleterious wtAA aaPos sfConsAA Pwt Psub message
 printf GOUT "$PTRout[0]\t$PTRout[1]\t$PTRout[2]\t$PTRout[3]\t$PTRout[4]\t$PTRout[13]\t$PTRout[14]\t$PTRout[15]\t$PTRout[20]\t";

#----------------------------process outputs------------------------------------
#Mutation Assesor
 my @MAout=open_file("MutationAssesor_output.txt");
 #legend
 #Func.Impact FI_score Func. region TS OG
 printf GOUT "$MAout[5]\t$MAout[6]\t$MAout[11]\t$MAout[12]\t$MAout[13]\t";

#PhD-SNP
 my @PhDout=open_file("PhD-SNP_output.txt");
 #legend
 #resn oAA nAA effect Rindex
 printf GOUT "$PhDout[4]\n";

#MAPP
 my @MAPPout=open_file("MAPP_output.txt");

#other
sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}

