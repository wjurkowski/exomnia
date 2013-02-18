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
my (@vpar,$fasta,$alignment,$tree);
my @inputs=open_file($ARGV[0]);

#output file
open (GOUT, ">SNPs_function_pipeline-out.txt") or die "Can not create general output file: $!";

#-----------------------create batch input for PPH and Panther---------------
#-----------------------run remaining tools-----------------------------------
open (MAOUT, "> MutationAssesor_output.txt") or die "Can not create MutAss output file: $!";
open (PDOUT, "> PhD-SNP_output.txt") or die "Can not create PhD-SNP output file: $!";
open (MAPPOUT, "> MAPP_output.txt") or die "Can not create MAPP output file: $!";
open (PTROUT, "> Panther_output.txt") or die "Can not create Panther output file: $!"; 
open (PPHOUT, "> Pph_output.txt") or die "Can not create PPH output file: $!";
open (PPH3D, "> Pph3D_output.txt") or die "Can not create PPH output file: $!";
print "start iterating input file--------------\n";
foreach my $lin (@inputs){
 @vpar=split(/\t/,$lin);
 $fasta=$vpar[0];
 my $uniprotAC=$vpar[1];
 my $uniprotID=$vpar[2];
 my $resn=$vpar[3];
 my $oAA=$vpar[4];
 my $nAA=$vpar[5];
 $nAA =~ s/\s+$//;
 $oAA =~ s/\s+$//;
 $alignment=$vpar[6];#in fasta format
 $tree=$vpar[7];
 my $waryjat="$oAA$resn$nAA";# gv e.g. G124A
 my @seq=open_file($fasta);

 #run Mutation Assesor
 print "run Mutation Assesor---------------------:\n";
 #batch mode by web api
 #required input: uniport ID and variation
 #http://mutationassessor.org/v1/?cm=var&p=EGFR_HUMAN&var=G719S
 my $URL="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
 my $response = $Agent->post($URL);
 sleep 5;
 die "$URL error: ", $response->status_line
 unless  $response->is_success;
 my $output="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
 my @out1 = get($output);
 my @MAo = split(/\t/,$out1[0]);
 printf MAOUT "$MAo[20]\t$MAo[21]\n";
 print "\tdone $waryjat on $uniprotAC\n";

 #run PhD-SNP
 print "run PhD-SNP------------------------------:\n";
 #python -O PhD-SNP.py -seq Test/1tthy.seq 21 K
 my $PhD="/usr/local/PhD-SNP2.0.6/";
 #my @seq_simple=grep(!/>/, @seq);
 `egrep -v '>' $fasta > sequence.txt`; 
 `python -O $PhD/PhD-SNP.py -seq sequence.txt $resn $nAA > outtmp`;
 my @out2=open_file("outtmp");
 my @PhDo=split(/ +/,$out2[10]);
 printf PDOUT "$PhDo[4]\t$PhDo[5]\n"; 
 print "\tdone $waryjat on $uniprotAC\n";

#run MAPP
print "run MAPP----------------------------------:\n";
#requires alignment in fasta format (ClustalW or ProbCons)
#requires tree build by ClustalW or Semphy
#output described in MAPP_readme.pdf
#physicochemical properties can be selected -s 1:3:4 (proterty 1+3+4), all on default 
#java -jar MAPP.jar -f LacI_Alignment.fa -t LacI.tree -o LacI_output.xls
#execute only if alignmnet and tree are present
 my $MAPP=$INST."MAPP";
 unless($alignment eq "-" and $tree eq "-"){
   `java -jar $MAPP/MAPP.jar -f $alignment -t $tree -o MAPP_tmp.txt`;
    my @out3=open_file("MAPP_tmp.txt");
    my @MAPPo=split(/ +/,$out3[0]);	
    printf MAPPOUT "$MAPPo[0]\n";
 }
 print "\tdone $waryjat on $uniprotAC\n";

#run Panther
print "run Panther-------------------------------:\n";
open (PTRIN, "> panther_csnpInput.txt") or die "Can not create Panther input file: $!"; 
printf PTRIN "$waryjat\|$uniprotAC\|$resn\|$oAA\;$nAA\n";
close (PTRIN);
#first classify protein sequence against PANTHER
 #pantherScore.pl -l <panther_hmm_library> -D B -V -i <fasta file> -o <output file> -n -T tmp/
 my $panther=$INST."Panther/csnpAnalysis1.02";
 my $tempdir="tmp";
 mkdir $tempdir;
 `$panther/pantherScore.pl -l $panther/PANTHER7.2 -D B -V -i $fasta -o panther_scores.out -n -T $tempdir`;
#classify SNPs
 #uprior.9comp contains some constants that could be modified
 #./snp_analysis.pl -l <panther_hmm_library> -c <score_outputfile_fromStepAbove> -s <csnp_input_file> -f <fasta_file> -b BLOSUM62 -V -p uprior.9comp  -o <output_file> -T tmp/ -a
 #If a SNP maps to multiple proteins, by default, the program will take the protein with the best HMm score.  Use the -a option if you want the results for all proteins
 #run analysis
 `$panther/snp_analysis.pl -l $panther/PANTHER7.2 -c panther_scores.out -s panther_csnpInput.txt -f $fasta -b $panther/BLOSUM62 -V -p $panther/uprior.9comp  -o panther_snpanalysis.out -T $tempdir`;
 `rm -rf $tempdir`;
#output format in README
 my @out5=open_file("panther_snpanalysis.out");
 for (my $i=1;$i<=$#out5;$i++){
  my @PTRo=split(/\t/,$out5[$i]);
  #output legend
  #snpId seqId subPSEC Pdeleterious wtAA aaPos sfConsAA Pwt Psub message
  if ($out5[$i] =~ m/.*SNP position in protein does not align to HMM/){
  	printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\tSNP outside HMM\n";
  	#print "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\tprobably non damaging\n";
  }
  elsif ($out5[$i] =~ m/.*invalid amino acid/){
  	printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\twrong input\n";
  	#print "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\t$PTRo[8]\n";
  }
  elsif ($out5[$i] =~ m/.*wild type amino acid is .*/){
  	printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\twrong input\n";
  	#print "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\t$PTRo[8]\n";
  }
  else{
  	printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\t$PTRo[5]\t$PTRo[6]\t$PTRo[11]\t$PTRo[12]\t$PTRo[13]\t$PTRo[14]\t$PTRo[19]\n";
  }
 }
print "\tdone $waryjat on $uniprotAC\n";

#run PolyPhen
print "run Polyphen------------------------------:\n";
open (PPHIN, "> pph_subs.input") or die "Can not create PPH input file: $!";
printf PPHIN "$uniprotAC\t$resn\t$oAA\t$nAA\n"; 
close (PPHIN);
#PolyPhen-2 analysis pipeline consists of three separate components,
 #each one executed by a dedicated Perl program:
 #
 #  * MapSNPs     (mapsnps.pl)   Genomic SNP annotation tool
 #  * PolyPhen-2  (run_pph.pl)   Protein variant annotation tool
 #  * PolyPhen-2  (run_weka.pl)  Probabilistic variant classifier
 #$PPH/bin/run_pph.pl subs.pph.input 1>pph.features 2>run_pph.log &
 my $PPH=$INST."PolyPhen/polyphen-2.2.2";
 `$PPH/bin/run_pph.pl pph_subs.input 1>pph.features 2>run_pph.log`;
 #$PPH/bin/run_weka.pl pph.features 1>pph.predictions
 `$PPH/bin/run_weka.pl pph.features 1>pph.predictions`;
 my @out4=open_file("pph.predictions");
 #my @out4n = grep(s/^\s*//g, @out4);
 #my (@out4n);
 #foreach my $el (@out4) {
 #  $el =~ s/^\s*(.*?)\s*$/$1/;
#print "$el\n";
 #  push (@out4n,$el);
 #}

# my $file="pph.predictions";
# open(INP, "< $file") or die "Can not open an input file: $!";
# my (@out4);
# while (<INP>) {
#  push(@out4, $_);
# }
 my @PPHo=split(/\t/,$out4[1]);
 s{^\s+|\s+$}{}g foreach @PPHo;

  #output legend
  #http://genetics.bwh.harvard.edu/pph2/dokuwiki/appendix_a
  #Uniprot_AC	tdbSNP	UniProtKB	pos	aa1	aa2	prediction	based_on	effect	pph2_class	pph2_prob	pph2_FPR	sensitivity	pph2_FDR	dScore	Score1	Score2	Nobs	Transv	CpG	MinDJnc	PfamHit
  printf PPHOUT "$PPHo[0]\t$PPHo[4]\t$PPHo[5]\t$PPHo[6]\t$PPHo[7]\t$PPHo[8]\t$PPHo[11]\t$PPHo[12]\t$PPHo[13]\t$PPHo[14]\t$PPHo[15]\t$PPHo[16]\t$PPHo[17]\t$PPHo[18]\t$PPHo[22]\t$PPHo[23]\t$PPHo[24]\t$PPHo[26]\t$PPHo[47]\t$PPHo[49]\t$PPHo[50]\t$PPHo[51]\n";
  #PDBID	PDB_pos	PDB_ch	ident	lenght	normSAS	SS	MapReg	dVol	dProp	B-fact	H-bonds	AveNHet	MinDHet	AveNInt	MinDInt	AveNSit	MinDSit
  printf PPH3D  "$PPHo[29]\t$PPHo[30]\t$PPHo[31]\t$PPHo[32]\t$PPHo[33]\t$PPHo[34]\t$PPHo[35]\t$PPHo[36]\t$PPHo[37]\t$PPHo[38]\t$PPHo[39]\t$PPHo[40]\t$PPHo[41]\t$PPHo[42]\t$PPHo[43]\t$PPHo[44]\t$PPHo[45]\t$PPHo[46]\n"; 
  print "\tdone $waryjat on $uniprotAC\n"; 

}# end of iteration
print "Iteration of input file done--------------\n";
close (PTROUT);
close (PPHOUT);
close (MAPPOUT);
close (PDOUT);
close (MAOUT);

#----------------------------process outputs------------------------------------
 my @out6=open_file("MutationAssesor_output.txt");
 my @out7=open_file("PhD-SNP_output.txt");
 my @out8=open_file("Panther_output.txt");
 my @out9=open_file("Pph_output.txt");
 my (@out10);
 unless($alignment eq "-" and $tree eq "-"){
  @out10=open_file("MAPP_output.txt");
 }
 for (my $i=0;$i<=$#inputs;$i++){
  #Panther
	my @PTRout=split(/\t+/,$out8[$i]);
	#legend
	#snpId seqId subPSEC Pdeleterious wtAA aaPos sfConsAA Pwt Psub message
	foreach my $val (@PTRout){
		printf GOUT "$val\t";
	}
  #Mutation Assesor
	my @MAout=split(/\t+/,$out6[$i]);
	#legend
	#Func.Impact FI_score Func. region TS OG
  	printf GOUT "$MAout[0]\t$MAout[1]\t";
  #PhD-SNP
	my @PhDout=split(/\s+/,$out7[$i]);
	#legend
	#resn oAA nAA effect Rindex
  	printf GOUT "$PhDout[0]\t$PhDout[1]\t";
  #Polyphen
	my @PPHout=split(/\t/,$out9[$i]);
	if($PPHout[8] =~ m/^$/){$PPHout[8]="effect unknown";}
	#legend
	#prediction	based_on	effect	pph2_prob pph2_FPR sensitivity
	printf GOUT "$PPHout[6]\t$PPHout[7]\t$PPHout[8]\t$PPHout[10]\t$PPHout[11]\t$PPHout[12]";
  #MAPP
 	unless($alignment eq "-" and $tree eq "-"){
	  my @MAPPout=split(/\s+/,$out10[$i]);
	  #legend
	}	
 printf GOUT "\n";#close current case
 }

#other
sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}


