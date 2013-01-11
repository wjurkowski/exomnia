#!/usr/bin/perl -w
use strict;
use warnings;
use LWP;
use LWP::Simple;
use LWP::UserAgent;
my $Agent = LWP::UserAgent->new;

if ($#ARGV != 1) {die "Program used with parameters [input file] [method] 

Available methods:
	0	all
	1	PhD-SNP
	2	Panther
	3	Polyphen
	4	SNAP
\n";}
#input file format: [fasta] [uniprot AC] [uniprot ID]\n"

#paths
my $INST="/home/wiktor/Komoda/Projects/Aktualne/Epilepsy/netmol_pipeline/multiplex/";
my $DB="$INST/benchmark/fasta/";

#input data
my (@vpar,$fasta,@AAsequence);
my $mALL=0;
my $mPhD=0;
my $mPTR=0;
my $mPPH=0;
my $mSNAP=0;
my @inputs=open_file($ARGV[0]);
@vpar=split(/\t/,$inputs[0]);
my @meth=split(/,/,$ARGV[1]);

foreach my $var (@meth){
	if($var == 0){$mALL=1;}
	if($var == 1){$mPhD=1;}
	if($var == 2){$mPTR=1;}
	if($var == 3){$mPPH=1;}	
	if($var == 4){$mSNAP=1;}
}


#process input sequence
$fasta="$DB/$vpar[0]";
my @seq=open_file($fasta);
foreach my $lin (@seq){
  unless($lin =~ m/>.*/){
    my @frag=split(//,$lin);
    push(@AAsequence,@frag);
  }
}
my $licznik = ($#AAsequence+1)*19;

#constants 
my(@AA);
$AA[0]="A";
$AA[1]="C";
$AA[2]="D";
$AA[3]="E";
$AA[4]="F";
$AA[5]="G";
$AA[6]="H";
$AA[7]="I";
$AA[8]="K";
$AA[9]="L";
$AA[10]="M";
$AA[11]="N";
$AA[12]="P";
$AA[13]="Q";
$AA[14]="R";
$AA[15]="S";
$AA[16]="T";
$AA[17]="V";
$AA[18]="W";
$AA[19]="Y";


#output files and input files for batch run
open (GOUT, ">SNPs_function_pipeline-out.txt") or die "Can not create general output file: $!";

if($mPhD == 1 or $mALL == 1){open (PDOUT, "> PhD-SNP_output.txt") or die "Can not create PhD-SNP output file: $!";}
if($mPTR == 1 or $mALL == 1){open (PTROUT, "> Panther_output.txt") or die "Can not create Panther output file: $!"; }
if($mPPH == 1 or $mALL == 1){
	open (PPHIN, "> pph_subs_input.txt") or die "Can not create PPH input file: $!";
	open (PPHOUT, "> PPH_output.txt") or die "Can not create PPH output file: $!";
	open (PPH3D, "> PPH3D_output.txt") or die "Can not create PPH output file: $!";
}
if($mSNAP == 1 or $mALL == 1){
	open (SNAPIN, "> snap_input.txt") or die "Can not create SNAP input file: $!";
	open (SNAPOUT, "> SNAP_output.txt") or die "Can not create SNAP output file: $!";
}

#-----------------------iterate over all positions and possible substitutions---------------
#


print "start iterating input file--------------\n";
my (@gener);
my $k=0;
foreach my $a (@AAsequence){
  $k++;
  for (my $i=0;$i<=19;$i++){
    unless($a eq $AA[$i]){
      my $uniprotAC=$vpar[1];
      my $uniprotID=$vpar[2];
      my $oAA=$a;
      my $resn=$k;
      my $nAA=$AA[$i];
      my $waryjat="$oAA$resn$nAA";# gv e.g. G124A
      my $zlepek="$uniprotAC\t$uniprotID\t$waryjat";
      push(@gener,$zlepek);

      if($mPhD == 1 or $mALL == 1){  #run PhD-SNP
	  print "run PhD-SNP------------------------------:\n";
	  #python -O PhD-SNP.py -seq Test/1tthy.seq 21 K
	  my $PhD="/usr/local/PhD-SNP2.0.6/";
	  #my @seq_simple=grep(!/>/, @seq);
	  `egrep -v '>' $fasta > sequence.txt`; 
	  `python -O $PhD/PhD-SNP.py -seq sequence.txt $resn $nAA > phd_tmp.out`;
	  my @PDout=open_file("phd_tmp.out");
	  my @PhDo=split(/ +/,$PDout[10]);
	  printf PDOUT "$PhDo[4]\t$PhDo[5]\n"; 
	  print "\tdone $waryjat on $uniprotAC\n";
	  `rm -f phd_tmp.out sequence.txt`;
      }

      if($mPTR == 1 or $mALL == 1){#run Panther
	  print "run Panther-------------------------------:\n";
	  open (PTRIN, "> panther_input.txt") or die "Can not create Panther input file: $!"; 
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
	  `$panther/snp_analysis.pl -l $panther/PANTHER7.2 -c panther_scores.out -s panther_input.txt -f $fasta -b $panther/BLOSUM62 -V -p $panther/uprior.9comp  -o panther_snpanalysis_tmp.out -T $tempdir`;
	  `rm -rf $tempdir`;
	  #output format in README
	  my @PTout=open_file("panther_snpanalysis_tmp.out");
	  for (my $i=1;$i<=$#PTout;$i++){
		my @PTRo=split(/\t/,$PTout[$i]);
		#output legend
		#snpId seqId subPSEC Pdeleterious wtAA aaPos sfConsAA Pwt Psub message
		if ($PTout[$i] =~ m/.*SNP position in protein does not align to HMM/){
			printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\tSNP outside HMM\n";
			#print "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\tprobably non damaging\n";
		}
		elsif ($PTout[$i] =~ m/.*invalid amino acid/){
			printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\twrong input\n";
			#print "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\t$PTRo[8]\n";
		}
		elsif ($PTout[$i] =~ m/.*wild type amino acid is .*/){
			printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\twrong input\n";
			#print "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\t$PTRo[8]\n";
		}
		else{
			printf PTROUT "$PTRo[0]\t$PTRo[1]\t$PTRo[2]\t$PTRo[3]\t$PTRo[4]\t$PTRo[5]\t$PTRo[6]\t$PTRo[11]\t$PTRo[12]\t$PTRo[13]\t$PTRo[14]\t$PTRo[19]\n";
		}
	  }
	  `rm -f panther_scores.out panther_snpanalysis_tmp.out panther_input.txt`;
	  print "\tdone $waryjat on $uniprotAC\n";
      }

#Polyphen-------------------------------------
      if($mPPH == 1 or $mALL == 1){#run PolyPhen	
	  #PolyPhen-2 analysis pipeline consists of three separate components,
	  #each one executed by a dedicated Perl program:
	  #
	  #  * MapSNPs     (mapsnps.pl)   Genomic SNP annotation tool
	  #  * PolyPhen-2  (run_pph.pl)   Protein variant annotation tool
	  #  * PolyPhen-2  (run_weka.pl)  Probabilistic variant classifier
	  printf PPHIN "$uniprotAC\t$resn\t$oAA\t$nAA\n"; 
      }

#SNAP-----------------------------------------
      if($mSNAP == 1 or $mALL == 1){#run SNAP
	  # snapinput format: AAposAA e.g C30Y
	  printf SNAPIN "$oAA$resn$nAA\n"; 
      }
 
    }# end of substitution iteration - synonymous skipped
  }# end of 20 aa iteration 
}# end of AA iteration
print "Iteration of sequence done--------------\n";

#close file handles and batch run
if($mPhD == 1 or $mALL == 1){close (PDOUT);}
if($mPTR == 1 or $mALL == 1){close (PTROUT);}
if($mPPH == 1 or $mALL == 1){
	print "run Polyphen------------------------------:\n";
	close (PPHIN);
	  my $PPH=$INST."PolyPhen/polyphen-2.2.2";
	  `$PPH/bin/run_pph.pl pph_subs_input.txt 1>pph_features.txt 2>pph.log`;
	  #$PPH/bin/run_weka.pl pph.features 1>pph.predictions
	  `$PPH/bin/run_weka.pl pph_features.txt 1>pph_predictions_tmp.out`;
	  my @PPHout=open_file("pph_predictions_tmp.out");
	  for(my $k=1;$k<=$#PPHout;$k++){
	  	my @PPHo=split(/\t/,$PPHout[$k]);
	  	s{^\s+|\s+$}{}g foreach @PPHo;#remove white spaces
	  	#output legend
	  	#http://genetics.bwh.harvard.edu/pph2/dokuwiki/appendix_a
	  	#Uniprot_AC	tdbSNP	UniProtKB	pos	aa1	aa2	prediction	based_on	effect	pph2_class	pph2_prob	pph2_FPR	sensitivity	pph2_FDR	dScore	Score1	Score2	Nobs	Transv	CpG	MinDJnc	PfamHit
	  	printf PPHOUT "$PPHo[0]\t$PPHo[4]\t$PPHo[5]\t$PPHo[6]\t$PPHo[7]\t$PPHo[8]\t$PPHo[11]\t$PPHo[12]\t$PPHo[13]\t$PPHo[14]\t$PPHo[15]\t$PPHo[16]\t$PPHo[17]\t$PPHo[18]\t$PPHo[22]\t$PPHo[23]\t$PPHo[24]\t$PPHo[26]\t$PPHo[47]\t$PPHo[49]\t$PPHo[50]\t$PPHo[51]\n";
	  	#PDBID	PDB_pos	PDB_ch	ident	lenght	normSAS	SS	MapReg	dVol	dProp	B-fact	H-bonds	AveNHet	MinDHet	AveNInt	MinDInt	AveNSit	MinDSit
	  	printf PPH3D  "$PPHo[29]\t$PPHo[30]\t$PPHo[31]\t$PPHo[32]\t$PPHo[33]\t$PPHo[34]\t$PPHo[35]\t$PPHo[36]\t$PPHo[37]\t$PPHo[38]\t$PPHo[39]\t$PPHo[40]\t$PPHo[41]\t$PPHo[42]\t$PPHo[43]\t$PPHo[44]\t$PPHo[45]\t$PPHo[46]\n"; 
	  	#`rm -f pph_features.txt pph_predictions_tmp.out pph.log pph_subs_input.txt`;
	  }
	close (PPHOUT);
	close (PPH3D);
}
if($mSNAP == 1 or $mALL == 1){
	print "run SNAP------------------------------:\n";
	close (SNAPIN);
	my $SNAP="/usr/bin";
	`$SNAP/snapfun -i $fasta -m snap_input.txt -o snap_tmp.out -a 20`;
	my @SNAPout=open_file("snap_tmp.out");
	for(my $k=12;$k<=$#SNAPout;$k++){
		my @SNAPo=split(/\t+/,$SNAPout[$k]);
		chomp @SNAPo;
		#output format
		#nsSNP   Prediction      Reliability Index       Expected Accuracy
		printf SNAPOUT "%s\t%s\t%s\t%s\n",$SNAPo[0],$SNAPo[1],$SNAPo[2],$SNAPo[3];
		#`rm -f snap_input.txt snap_tmp.out`;
	}
	close (SNAPOUT);
}

#----------------------------process outputs------------------------------------
my (@out5,@out7,@out8,@out9);
if($mPhD == 1 or $mALL == 1){@out7=open_file("PhD-SNP_output.txt");}
if($mPTR == 1 or $mALL == 1){@out8=open_file("Panther_output.txt");}
if($mPPH == 1 or $mALL == 1){@out9=open_file("PPH_output.txt");}
if($mSNAP == 1 or $mALL == 1){@out5=open_file("SNAP_output.txt");}

 
 for (my $i=0;$i<=$licznik-1;$i++){
	printf GOUT "$gener[$i]\t";
	if($mPTR == 1 or $mALL == 1){
		#Panther
		my @PTRout=split(/\t+/,$out8[$i]);
		#legend
		#snpId seqId subPSEC Pdeleterious wtAA aaPos sfConsAA Pwt Psub message
		for (my $val=3;$val <= $#PTRout;$val++){
			printf GOUT "$PTRout[$val]\t";
		}
	}	
	if($mPhD == 1 or $mALL == 1){
		#PhD-SNP
		my @PhDout=split(/\s+/,$out7[$i]);
		#legend
		#resn oAA nAA effect Rindex
		printf GOUT "$PhDout[0]\t$PhDout[1]\t";
	}
	if($mPPH == 1 or $mALL == 1){	
		#Polyphen
		my @PPHout=split(/\t/,$out9[$i]);
		if($PPHout[8] =~ m/^$/){$PPHout[8]="effect unknown";}
		#legend
		#prediction	based_on	effect	pph2_prob pph2_FPR sensitivity
		printf GOUT "$PPHout[6]\t$PPHout[7]\t$PPHout[8]\t$PPHout[10]\t$PPHout[11]\t$PPHout[12]\t";
	}
	if($mSNAP == 1 or $mALL == 1){	
		#SNAP
		my @SNAPout=split(/\t/,$out5[$i]);
		printf GOUT "%s\t%s\t%s\t",$SNAPout[1],$SNAPout[2],$SNAPout[3];

	}
	
 printf GOUT "\n";#close current case
 }#end if GV iteration

#other
sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}


