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
	1	MutationAssesor
	2	PhD-SNP
	3	MAPP
	4	Panther
	5	Polyphen
	6	SNAP
\n";}
#input file format: [fasta] [uniprot AC] [resn] [oAA] [nAA] [alignment] [tree]\n"

#paths
my $INST="/home/wiktor/Komoda/Projects/Aktualne/Epilepsy/netmol_pipeline/multiplex/";
my $DB="$INST/benchmark/fasta/";

#input data
my (@vpar,$fasta,$alignment,$tree);
my $mALL=0;
my $mMA=0;
my $mPhD=0;
my $mMAPP=0;
my $mPTR=0;
my $mPPH=0;
my $mSNAP=0;
my $mNET=0;
my @inputs=open_file($ARGV[0]);
my @meth=split(/,/,$ARGV[1]);
foreach my $var (@meth){
	if($var == 0){$mALL=1;}
	if($var == 1){$mMA=1;}
	if($var == 2){$mPhD=1;}
	if($var == 3){$mMAPP=1;}
	if($var == 4){$mPTR=1;}
	if($var == 5){$mPPH=1;}	
	if($var == 6){$mSNAP=1;}
	if($var == 7){$mNET=1;}	
}

#output file
open (GOUT, ">SNPs_function_pipeline-out.txt") or die "Can not create general output file: $!";

if($mMA == 1 or $mALL == 1){open (MAOUT, "> MutationAssesor_output.txt") or die "Can not create MutAss output file: $!";}
if($mPhD == 1 or $mALL == 1){open (PDOUT, "> PhD-SNP_output.txt") or die "Can not create PhD-SNP output file: $!";}
if($mMAPP == 1 or $mALL == 1){open (MAPPOUT, "> MAPP_output.txt") or die "Can not create MAPP output file: $!";}
if($mPTR == 1 or $mALL == 1){open (PTROUT, "> Panther_output.txt") or die "Can not create Panther output file: $!"; }
if($mPPH == 1 or $mALL == 1){
	open (PPHOUT, "> PPH_output.txt") or die "Can not create PPH output file: $!";
	open (PPH3D, "> PPH3D_output.txt") or die "Can not create PPH output file: $!";
}
if($mSNAP == 1 or $mALL == 1){open (SNAPOUT, "> SNAP_output.txt") or die "Can not create SNAP output file: $!";}

#-----------------------run SNP effect predictors---------------
#
my (@gener,$l,$runfasta);
my $prot="blaszka";

print "start iterating input file--------------\n";
for ($l=0;$l<=$#inputs;$l++){
 @vpar=split(/\t/,$inputs[$l]);
 $fasta="$DB/$vpar[0]";
 my $uniprotAC=$vpar[1];
 my $uniprotID=$vpar[2];
 my $oAA=$vpar[3];
 my $resn=$vpar[4];
 my $nAA=$vpar[5];
 $nAA =~ s/\s+$//;
 $oAA =~ s/\s+$//;
 $alignment=$vpar[6];#in fasta format
 $tree=$vpar[7];
 my $waryjat="$oAA$resn$nAA";# gv e.g. G124A
 my @seq=open_file($fasta);
 my $zlepek="$uniprotAC\t$uniprotID\t$waryjat";
 push(@gener,$zlepek);

 if($mMA == 1 or $mALL == 1){
	#run Mutation Assesor
	print "run Mutation Assesor---------------------:\n";
	#batch mode by web api
	#required input: uniport ID and variation
	#http://mutationassessor.org/v1/?cm=var&p=EGFR_HUMAN&var=G719S
	my $output="http://mutationassessor.org/v1/?cm=var&p=$uniprotID&var=$waryjat&frm=txt";
	my $str = get($output);
	$str =~ s/\n/ /g;
	my @MAo = split(/\t/,$str);
	unless($MAo[18] =~ /http/){
	  $MAo[20]="none";
	  $MAo[21]=-1;
	}
	printf MAOUT "$MAo[20]\t$MAo[21]\n";
	print "\tdone $waryjat on $uniprotAC\n";
}

 if($mPhD == 1 or $mALL == 1){
	#run PhD-SNP
	print "run PhD-SNP------------------------------:\n";
	#python -O PhD-SNP.py -seq Test/1tthy.seq 21 K
	my $PhD="/usr/local/PhD-SNP2.0.6/";
	#my @seq_simple=grep(!/>/, @seq);
	`egrep -v '>' $fasta > sequence.txt`; 
	`python -O $PhD/PhD-SNP.py -seq sequence.txt $resn $nAA > phd_tmp.out`;
	my @PDout=open_file("phd_tmp.out");
	my @PhDo=split(/ +/,$PDout[10]);
	if($PhDo[4] =~ m/^$/){$PhDo[4]="none";}
	if($PhDo[5] =~ m/^$/){$PhDo[5]="-1";}
	printf PDOUT "$PhDo[4]\t$PhDo[5]\n"; 
	print "\tdone $waryjat on $uniprotAC\n";
	`rm -f phd_tmp.out sequence.txt`;
 }

 if($mMAPP == 1 or $mALL == 1){
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
		`java -jar $MAPP/MAPP.jar -f $alignment -t $tree -o mapp_tmp.out`;
		my @MPout=open_file("mapp_tmp.out");
		my @MAPPo=split(/ +/,$MPout[0]);	
		printf MAPPOUT "$MAPPo[0]\n";
	}
	print "\tdone $waryjat on $uniprotAC\n";
	`rm -f mapp_tmp.out`;
 } 

 if($mPTR == 1 or $mALL == 1){
	#run Panther
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

 if($mPPH == 1 or $mALL == 1){
	#run PolyPhen
	my @PPHo;
	print "run Polyphen------------------------------:\n";
	open (PPHIN, "> pph_subs_input.txt") or die "Can not create PPH input file: $!";
	printf PPHIN "$uniprotAC\t$resn\t$oAA\t$nAA\n"; 
	close (PPHIN);
	#PolyPhen-2 analysis pipeline consists of three separate components,
	#each one executed by a dedicated Perl program:
	#
	#  * MapSNPs     (mapsnps.pl)   Genomic SNP annotation tool
	#  * PolyPhen-2  (run_pph.pl)   Protein variant annotation tool
	#  * PolyPhen-2  (run_weka.pl)  Probabilistic variant classifier
	#$PPH/bin/run_pph.pl subs.pph.input 1>pph.features 2>pph.log &
	my $PPH=$INST."PolyPhen/polyphen-2.2.2";
	`$PPH/bin/run_pph.pl pph_subs_input.txt 1>pph_features.txt 2>pph.log`;
	#$PPH/bin/run_weka.pl pph.features 1>pph.predictions
	`$PPH/bin/run_weka.pl pph_features.txt 1>pph_predictions_tmp.out`;
	my @PPout=open_file("pph_predictions_tmp.out");
	@PPHo=split(/\t/,$PPout[1]) if (exists $PPout[1]);
	s{^\s+|\s+$}{}g foreach @PPHo;#remove white spaces
	if(exists $PPHo[0]){
	  if($PPHo[6] =~ m/^$/){$PPHo[6]="none";}
	  if($PPHo[7] =~ m/^$/){$PPHo[7]="none";}
	  if($PPHo[8] =~ m/^$/){$PPHo[8]="effect unknown";}
	  if($PPHo[10] =~ m/^$/){$PPHo[10]=-1;}
	  if($PPHo[11] =~ m/^$/){$PPHo[11]=-1;}
	  if($PPHo[12] =~ m/^$/){$PPHo[12]=-1;}
	}
	else{
	  $PPHo[6]="no res";
	  $PPHo[7]="no res";
	  $PPHo[8]="no res";
	  $PPHo[10]="no res";
	  $PPHo[11]="no res";
	  $PPHo[12]="no res"; 
	}
	if(exists $PPHo[0]){
	  #output legend
	  #http://genetics.bwh.harvard.edu/pph2/dokuwiki/appendix_a
	  #Uniprot_AC	tdbSNP	UniProtKB	pos	aa1	aa2	prediction	based_on	effect	pph2_class	pph2_prob	pph2_FPR	sensitivity	pph2_FDR	dScore	Score1	Score2	Nobs	Transv	CpG	MinDJnc	PfamHit
	  printf PPHOUT "$PPHo[0]\t$PPHo[4]\t$PPHo[5]\t$PPHo[6]\t$PPHo[7]\t$PPHo[8]\t$PPHo[11]\t$PPHo[12]\t$PPHo[13]\t$PPHo[14]\t$PPHo[15]\t$PPHo[16]\t$PPHo[17]\t$PPHo[18]\t$PPHo[22]\t$PPHo[23]\t$PPHo[24]\t$PPHo[26]\t$PPHo[47]\t$PPHo[49]\t$PPHo[50]\t$PPHo[51]\n";
	  #PDBID	PDB_pos	PDB_ch	ident	lenght	normSAS	SS	MapReg	dVol	dProp	B-fact	H-bonds	AveNHet	MinDHet	AveNInt	MinDInt	AveNSit	MinDSit
	  printf PPH3D  "$PPHo[29]\t$PPHo[30]\t$PPHo[31]\t$PPHo[32]\t$PPHo[33]\t$PPHo[34]\t$PPHo[35]\t$PPHo[36]\t$PPHo[37]\t$PPHo[38]\t$PPHo[39]\t$PPHo[40]\t$PPHo[41]\t$PPHo[42]\t$PPHo[43]\t$PPHo[44]\t$PPHo[45]\t$PPHo[46]\n"; 
	  print "\tdone $waryjat on $uniprotAC\n"; 
	  #`rm -f pph_features.txt pph_predictions_tmp.out pph.log pph_subs_input.txt`;
	}
	else{
	  printf PPHOUT "ERROR: no results\n";
	  printf PPH3D "ERROR: no results\n"
	}
 }

 if($mSNAP == 1 or $mALL == 1){
	#run SNAP
	my $SNAP="/usr/local/bin";
	mkdir "SNAP_tmpdir";

	if($uniprotAC eq $prot){
	  # snapinput format: AAposAA e.g C30Y
	  printf SNAPIN "$oAA$resn$nAA\n"; 
	}
	else{
	  if($prot eq "blaszka"){
		open (SNAPIN, "> snap_input.txt") or die "Can not create SNAP input file: $!";	
	  	printf SNAPIN "$oAA$resn$nAA\n"; 
	  }
	  unless($prot eq "blaszka"){
	  	close (SNAPIN);
		print "run SNAP------------------------------:\n";
		`$SNAP/slap -i $runfasta -m snap_input.txt -o snap_tmp.out -w SNAP_tmpdir -cpus 20`;
		my @SNAPout=open_file("snap_tmp.out");
		for(my $k=12;$k<=$#SNAPout;$k++){
	  		my @SNAPo=split(/\s+/,$SNAPout[$k]);
	  		chomp @SNAPo;
	  		#output format
	  		#nsSNP   Prediction      Reliability Index       Expected Accuracy
	  		printf SNAPOUT "%s\t%s\t%s\t%s\n",$SNAPo[0],$SNAPo[1],$SNAPo[2],$SNAPo[3];
	  		#`rm -f snap_input.txt snap_tmp.out`;
	  	}
	  	open (SNAPIN, "> snap_input.txt") or die "Can not create SNAP input file: $!";
	  	printf SNAPIN "$oAA$resn$nAA\n"; 
	  }
	}
	$prot = $uniprotAC;
	$runfasta = $fasta;
	if($l == $#inputs){
	  	close (SNAPIN);
		print "run SNAP------------------------------:\n";
		`$SNAP/slap -i $fasta -m snap_input.txt -o snap_tmp.out -w SNAP_tmpdir -cpus 20`;
		my @SNAPout=open_file("snap_tmp.out");
		for(my $k=12;$k<=$#SNAPout;$k++){
	  		my @SNAPo=split(/\s+/,$SNAPout[$k]);
	  		chomp @SNAPo;
		  	#output format
		  	#nsSNP   Prediction      Reliability Index       Expected Accuracy
	  		printf SNAPOUT "%s\t%s\t%s\t%s\n",$SNAPo[0],$SNAPo[1],$SNAPo[2],$SNAPo[3];
		  	#`rm -f snap_input.txt snap_tmp.out`;
	  	}
	}
	#about SNAP
	#cleaning tmp files
	#`rm -f SNAP_tmpdir/*`;
	#rmdir ("SNAP_tmpdir") || die ("error in deleting directory: $?");
 }
	
 if($mNET == 1 or $mALL == 1){	
	#grab network scores
	print "grab network---------------------------:\n";
	#no need for reformating now
 }

}# end of iteration
print "Iteration of input file done--------------\n";


#close file handles
if($mMA == 1 or $mALL == 1){close (MAOUT);}
if($mPhD == 1 or $mALL == 1){close (PDOUT);}
if($mMAPP == 1 or $mALL == 1){close (MAPPOUT);}
if($mPTR == 1 or $mALL == 1){close (PTROUT);}
if($mPPH == 1 or $mALL == 1){
	close (PPHOUT);
	close (PPH3D);
}
if($mSNAP == 1 or $mALL == 1){close (SNAPOUT);}

#----------------------------process outputs------------------------------------
my (@out4,@out5,@out6,@out7,@out8,@out9,@out10);
if($mMA == 1 or $mALL == 1){@out6=open_file("MutationAssesor_output.txt");}
if($mPhD == 1 or $mALL == 1){@out7=open_file("PhD-SNP_output.txt");}
if($mMAPP == 1 or $mALL == 1){@out10=open_file("MAPP_output.txt");}
if($mPTR == 1 or $mALL == 1){@out8=open_file("Panther_output.txt");}
if($mPPH == 1 or $mALL == 1){@out9=open_file("PPH_output.txt");}
if($mSNAP == 1 or $mALL == 1){@out5=open_file("SNAP_output.txt");}
if($mNET == 1 or $mALL == 1){@out4=open_file("network_scores.txt");}

 for (my $i=0;$i<=$#inputs;$i++){
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
	if($mMA == 1 or $mALL == 1){	
		#Mutation Assesor
		my @MAout=split(/\t+/,$out6[$i]);
		#legend
		#Func.Impact FI_score Func. region TS OG
		printf GOUT "$MAout[0]\t$MAout[1]\t";
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
		unless($out9[0] =~ /ERROR/){
		  my @PPHout=split(/\t/,$out9[$i]);
		  #legend
		  #prediction	based_on	effect	pph2_prob pph2_FPR sensitivity
		  printf GOUT "$PPHout[6]\t$PPHout[7]\t$PPHout[8]\t$PPHout[10]\t$PPHout[11]\t$PPHout[12]\t";
		}
	}
	if($mMAPP == 1 or $mALL == 1){	
		#MAPP
		my @MAPPout=split(/\s+/,$out10[$i]);
		#legend
	}
	if($mSNAP == 1 or $mALL == 1){	
		#SNAP
		my @SNAPout=split(/\t/,$out5[$i]);
		printf GOUT "%s\t%s\t%s\t",$SNAPout[1],$SNAPout[2],$SNAPout[3];
	}
	if($mNET == 1 or $mALL == 1){
		#net
		my @Netout=split(/\s+/,$out4[$i]);
		foreach my $val (@Netout){
			printf GOUT "$val\t";
		}
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


