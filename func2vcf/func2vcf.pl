# Program translating func file (CCG Cologne) to VCF4.1 format
# # Copyright:: Copyright 2012-2013 University of Luxembourg/Luxembourg Centre for Systems Biomedicine
# Authors::   Wiktor Jurkowski <wiktor.jurkowski@uni.lu>
# WARNING: called variants with zero frequency and noncalled positions are not sorted out

#!/usr/bin/perl -w
use strict;
use warnings;
#call external modules
use Math::Round qw(nearest);
use Time::Piece;
use Getopt::Std;

#control run options
my (%opt);
getopts('htvf:o:',\%opt);
die &usage() if (tutej_any($opt{f},$opt{t},$opt{v},$opt{o}));
&usage() if $opt{h};

#open output file
open(GOUT, ">$opt{o}.vcf") or die "Can not open an output file: $!";#output file
my $stats=$opt{o}.".QAsummary.tsv";
my $stats2=$opt{o}.".QAsummary_scoresSNP.tsv";
my $stats3=$opt{o}.".QAsummary_scoresINDEL.tsv";
my $stats_all=$opt{o}.".ALLstats.tsv";
my $stats_gatk=$opt{o}.".GATKstats.tsv";
my $stats_gatk2=$opt{o}.".GATKstats-solo.tsv";
my $stats_mpileup=$opt{o}.".MPILEUPstats.tsv";
my $stats_mpileup2=$opt{o}.".MPILEUPstats-solo.tsv";
my $stats_dindel=$opt{o}.".DINDELstats.tsv";
my $stats_dindel2=$opt{o}.".DINDELstats-solo.tsv";
open (GSTAT,">$stats") or die "Can not open an output file: $!";
open (GSNP,">$stats2") or die "Can not open an output file: $!";
open (GIND,">$stats3") or die "Can not open an output file: $!";
open (ALL,">$stats_all") or die "Can not open an output file: $!";
open (GATK,">$stats_gatk") or die "Can not open an output file: $!";
open (GATK2,">$stats_gatk2") or die "Can not open an output file: $!";
open (MPILEUP,">$stats_mpileup") or die "Can not open an output file: $!";
open (MPILEUP2,">$stats_mpileup2") or die "Can not open an output file: $!";
open (DINDEL,">$stats_dindel") or die "Can not open an output file: $!";
open (DINDEL2,">$stats_dindel2") or die "Can not open an output file: $!";
printf GSTAT "Sample\tTotal\tGATK\tMPILEUP\tDINDEL\tMRSSVSP\tGATKs\tMPILEUPs\tDINDELs\tMRSSVSPs\n";
printf GSNP "Sample\tpass0\tpass1\tpass2\tpass3\tpass4\tpass5\tpass6\tQDpass\tQDfail\tMQpass\tMQfail\tFSpass\tFSfail\tHSpass\tHSfail\tMQRSpass\tMQRSfail\tRPRSpass\tRPRSfail\tVQSLODpass\tVQSLODfail\n";
printf GIND "Sample\tpass0\tpass1\tpass2\tpass3\tpass4\tQDpass\tQDfail\tFSpass\tFSfail\tRPRSpass\tRPRSfail\tICpass\tICfail\n";
printf ALL "Sample\tType\tGATK\tMPILEUP\tDINDEL\tCoverage\tRefF\tAltF\tGQ\tQD\tMQ\tFS\tHaplotypeScore\tMQRankSum\tReadPosRankSum\tBaseQRankSum\tVQSLOD\tIC\tMMQ\tFQ\tPV0\tPV1\tPV2\tPV3\tPL0\tPL1\tPL2\tSP\tDGQ\tNF\tNR\tNRS\tNFS\tHP\n";
printf GATK "Sample\tCaller\tType\tCoverage\tRefF\tAltF\tGQ\tQD\tMQ\tFS\tHaplotypeScore\tMQRankSum\tReadPosRankSum\tBaseQRankSum\tVQSLOD\tIC\n";
printf GATK2 "Sample\tCaller\tType\tCoverage\tRefF\tAltF\tGQ\tQD\tMQ\tFS\tHaplotypeScore\tMQRankSum\tReadPosRankSum\tBaseQRankSum\tVQSLOD\tIC\n";
printf MPILEUP "Sample\tCaller\tType\tCoverage\tRefF\tAltF\tMQ\tFQ\tPV0\tPV1\tPV2\tPV3\tPL0\tPL1\tPL2\tSP\n";
printf MPILEUP2 "Sample\tCaller\tType\tCoverage\tRefF\tAltF\tMQ\tFQ\tPV0\tPV1\tPV2\tPV3\tPL0\tPL1\tPL2\tSP\n";
printf DINDEL "Sample\tCaller\tType\tCoverage\tRefF\tAltF\tGQ\tNF\tNR\tNRS\tNFS\tHP\n";
printf DINDEL2 "Sample\tCaller\tType\tCoverage\tRefF\tAltF\tGQ\tNF\tNR\tNRS\tNFS\tHP\n";


#define time
my $date = Time::Piece->new->strftime('%m/%d/%Y');

#print header
print GOUT "##fileformat=VCFv4.1\n";
print GOUT "##fileDate=$date\n";
print GOUT "##source=CCG func files\n";
print GOUT "##reference=NCBI37\n";
print GOUT "##phasing=partial\n";
print GOUT "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
print GOUT "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
print GOUT "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes\">\n";
print GOUT "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
print GOUT "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n";
print GOUT "##INFO=<ID=HGNC,Number=.,Type=String,Description=\"HGNC gene symbol\">\n";
print GOUT "##INFO=<ID=CDNA,Number=.,Type=String,Description=\"cDNA variation\">\n";
print GOUT "##INFO=<ID=VARAA,Number=.,Type=String,Description=\"AA variation\">\n";
print GOUT "##INFO=<ID=TYPE,Number=.,Type=String,Description=\"Variation type: SNP, CNV, INDEL etc.\">\n";
print GOUT "##INFO=<ID=CONSQ,Number=.,Type=String,Description=\"Consequence: synonymous, non-synonymous, splicing etc.\">\n";
print GOUT "##INFO=<ID=ET,Number=.,Type=String,Description=\"ENSEMBL transcript ID\">\n";
print GOUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print GOUT "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotyping Quality\">\n";
print GOUT "##FORMAT=<ID=VQ,Number=1,Type=Float,Description=\"VQSLOD\">\n";
print GOUT "##FORMAT=<ID=CP,Number=1,Type=Integer,Description=\"GATK criteria passed\">\n";
print GOUT "##FORMAT=<ID=BF,Number=1,Type=Float,Description=\"Alt base frequency\">\n";
print GOUT "##FORMAT=<ID=NC,Number=1,Type=Integer,Description=\"Number of callers\">\n";
print GOUT "##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Read Depth\">\n";
print GOUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";#header of fixed1

#read list of func files
open(DANE, "< $opt{f}") or die "Can not open an input file: $!";#reads in the file
my @flist=<DANE>;
close (DANE);
chomp @flist;

#define global variables
my (%AC1,%AC2,%genotyping,%fixed1,$format,%an1,%an2,%an3,%an4,%an5,%an6,%alts,@sex);
#read func file
my $NS=0;
my $males=0;
foreach my $funcf(@flist){
  ###########################################################
  ##initial for each data set
  $NS++;
  open(FUNC, "< $funcf") or die "Can not open an input file: $!";#reads in the func file
  my $funcn=$funcf;
  $funcn=~s/\.func//;
  $funcn=$funcn.".clean.func";
  open (FUNCN,">$funcn") or die "Can not open an output file: $!";#open clean func file
  $funcf=~s/.*\///;
  $funcf=~s/\..*.func//;
  my $errorf=$funcf.".except.tsv";
  open (FUNCE,">$errorf") or die "Can not open an output file: $!";#open error file
  $funcf=~s/AID\d+_//;
  $funcf=~s/_(.*)$//;
  print GOUT "\t$funcf"; #data set ID
  my (%genotype,@totcdna,@totprot,@totenst,@totconsq);
  my $known="NA";
  my $is_male = my $gatk = my $mpileup = my $dindel = my $mrssvsp = my $gsolo = my $msolo = my $dsolo = my $mrsolo = 0;
  my $gatk0 = my $gatk1 = my $gatk2 = my $gatk3 = my $gatk4 = my $gatk5 = my $gatk6 = my $igatk0 = my $igatk1 = my $igatk2 = my $igatk3 = my $igatk4 = 0;
  my $qdpass = my $qdfail = my $mqpass = my $mqfail = my $fspass = my $fsfail = my $hspass = my $hsfail = 0;
  my $mqrspass = my $mqrsfail = my $rprspass = my $rprsfail = my $vqslodpass = my $vqslodfail = 0;
  my $iqdpass = my $iqdfail =  my $ifspass = my $ifsfail = my $irprspass = my $irprsfail = my $icpass = my $icfail = my $ivqslodpass = my $ivqslodfail = 0;
  ###########################################################
  ##check sex
  while(<FUNC>){
    if(/^Y/){
      $is_male=1;
      last;
    }
  }
  seek(FUNC,0,0);

  ##read in exons for sample in question
  #my (@exons);
  #if($opt{t}){
  #  my @exonf = grep /$funcf/, @elist;
  #  open (ECOV,"<$exonf[0]") or die "Can not open an input file: $!";#open exon file
  #  @exons=<ECOV>;
  #  chomp @exons;
  #  close (ECOV);
  #}

  while(<FUNC>){#for variant in each data set
   ###########################################################
   ##reads, initial quality check
   my $variant = $_;
   chomp $variant;
   unless($variant =~ /^#/){
    my @data=split("\t",$variant);
    #test variation type
    my $type=$data[5];
    my $target=$data[10];
    my $chrom=$data[0];

    ##QC step 1
    #skip ROH variants
    if($type eq "ROH" or $type eq "CNV"){
	print FUNCE "ROH_CNV\t$variant\n";
	next;
    }
    ##QC step 2
    #QA: move lines with uncomplete data to error file
    #skip lines with missing or ambigous genotype info
    if($data[13] eq "" or $data[13] eq "?"){
	print FUNCE "NO_GENOTYPE\t$variant\n";
	next;
    }
    ##QC step 3
    #skip variant outside of target regions (+/- 100bp)
    if($opt{t}){
      if($target){
	if($target > 100){
	  next;
	}
      }
      else{
	next unless $chrom eq "MT";
      }
    }

    ###########################################################
    ##read in data
    my $altF = my $refF = "NA";
    my $pos=$data[1];
    my $ref=$data[3];
    my $alt=$data[4];
    my $qual=$data[9];
    $qual="." if $qual eq "";
    my $filter=".";#default no filter info 
    $altF=$data[14];
    $refF=$data[15];
    my $id=".";#missing on default
    $id=$data[16] if ($data[16]);#dbSNP or other ID
    my $cover=$data[8];
    my $ensg=$data[21];
    my $enst=$data[22];
    my $ensp=$data[23];
    my $biot=$data[25];
    my $hgnc=$data[26];
    my $region=$data[32];
    my $consq=$data[33];
    $consq=~s/;/,/g;
    my $cdna=$data[34];
    my $prot=$data[35];

    ##########################################################
    ##interpret caller info: 
    $data[20] =~ s/GATK=/GATK:/;#remove error of incorectly formatted GATK caller info
    my @callers=split(']',$data[20]);
    my $cn=@callers; 
    my $isgatk = my $isdindel = my $ismpileup = my $ismrssvsp = 0;
    my (@gloc,@mloc,@dloc);
    foreach my $c(@callers){
	$c =~ s/;\[//;#remove caller info brackets 
	$c =~ s/\[//;
	my (%stats);
	my @caller=split(":",$c);
	my $calln=shift(@caller);
	my $rest=join(':',@caller);
	my @callerinfo=split(";",$rest);
	foreach my $field(@callerinfo){
	  my @pair=split("=",$field);
         $stats{$pair[0]}=$pair[1];
	}
	my $locorg=$stats{"LOC_ORG"} if $stats{"LOC_ORG"};
	#go through callers	
 	if($calln eq "GATK" or $calln eq "GATK_SANGER"){
	  @gloc=split(",",$locorg);# original locus info
	  $isgatk = 1;
	}
	elsif($calln eq "MPILEUP"){
	  @mloc=split(",",$locorg);# original locus info
	  $ismpileup = 1;
	}
	elsif($calln eq "DINDEL" or $calln eq "DINDEL_SANGER"){
	  @dloc=split(",",$locorg);# original locus info
	  $isdindel = 1;
	}
	elsif($calln eq "MRSSVSP"){
	   $ismrssvsp = 1;
	}
    }
    #oryginal allele info 
    #is taken from last caller on the list: should be same in all
    my ($opos,$oref,$oalt,$nalt);
    my $noloc=0;
    if($isgatk == 1){
      if(@gloc){
	$opos=shift(@gloc);#position
	$oref=shift(@gloc);#reference allele
	$nalt=@gloc;
	$oalt=join(",",@gloc);#alternative allele
	$oalt=~s/\]//;#remove closing bracket if at the end of caller field
      }
      else{
	$noloc=1;
      }
    }
    else{#not gatk
      if($ismpileup == 1){
	if(@mloc){
	  $opos=shift(@mloc);#position
	  $oref=shift(@mloc);#reference allele
	  $nalt=@mloc;
	  $oalt=join(",",@mloc);#alternative allele
	  $oalt=~s/\]//;#remove closing bracket if at the end of caller field
	}
	else{
	  $noloc=1;
	}
      }
      else{#not mpileup
	if($isdindel == 1){
	 if(@dloc){
	   my @dloc2 = grep ! /<DEL>/, @dloc;
	   $opos=shift(@dloc2);#position
	   $oref=shift(@dloc2);#reference allele
	   $nalt=@dloc2;
	   $oalt=join(",",@dloc2);#alternative allele
	   $oalt=~s/\]//;#remove closing bracket if at the end of caller field
	 }
	 else{
	   $noloc=1;
	 }
	}
	else{#not dindel
	 if($ismrssvsp == 1){
	  $opos=$pos-1;
	  $nalt=1;
	  if($ref=~/\(/){
		$ref=~s/\(//;
		$ref=~s/\)/;/;
		my @hind=split(";",$ref);
		$ref=~s/$hind[0]\;\d+/$hind[0] x $hind[1]/e;
	  }
	  $oref=$ref;
	  my @mrloc=split(",",$alt);
	  $nalt=@mrloc;
	  my @corr;
	  foreach my $v(@mrloc){
	    if($v=~/\(/){
		$v=~s/\(//;
		$v=~s/\)/;/;
		my @hind=split(";",$v);
		$v=~s/$hind[0]\;\d+/$hind[0] x $hind[1]/e;
		push(@corr,$v);
	    }
	  }
	  $oalt=join(",",@corr);
	 }
	 else{#not mrssvsp
	   $noloc=1;
	   print "przeca $variant\n";
	 }
	}
      }
    }#end of original caller info check
    my $var=$chrom."-".$opos.$oalt; #each variation is a key
  
    ###########################################################
    ##QC step 4
    #skip lines without LOC_ORG (except MRSSVSP that is consistently lacking it)
    unless($noloc == 0){
	print FUNCE "NO_LOC_ORG\t$variant\n";
	next;
    }
    ##skip correct redundant lines annotated with non-coding transcript ENSEMBLT:
    unless($biot eq "protein_coding"){
	print FUNCE "NON PROTEIN CODING TRANSCRIPT\t$variant\n";
	next;
    }
    ##QC step 5  
    #correct incorrect delimiter in genotype field
    $data[13]=~s/:/\//;
    #remove variants with mismatching alt alleles # and genotype
    if($nalt > 1 and ($data[13] eq "0|1" or $data[13] eq "0/1" or $data[13] eq "1|0" or $data[13] eq "1|1" or $data[13] eq "1/1")){
	print FUNCE "WRONG_GENOTYPE_OR_ALT: Number of alt alleles higher than expected for given genotype\t$variant\n";
	next;
    }
    elsif($nalt == 1 and ($data[13] eq "1|2" or $data[13] eq "1/2" or $data[13] eq "2|1" or $data[13] eq "2|2" or $data[13] eq "2/2" or $data[13] eq "0|2" or $data[13] eq "0/2" or $data[13] eq "2|0" )){
	print FUNCE "WRONG_GENOTYPE_OR_ALT: Number of alt alleles lower than expected for given genotype\t$variant\n";
	next;
    }
    if($is_male == 1 and ($chrom eq "X" or $chrom eq "Y")){
	if($nalt > 1){
		print FUNCE "WRONG_GENOTYPE_OR_ALT: Number of alt alleles higher than expected for given genotype\t$variant\n";
		next;
	}
    }
    if($data[12] eq 2 and $data[13] eq "1/1"){
	print FUNCE "REDUNDANT ENTRY; Allele ID =2 but given genotype has only one allele type\t$variant\n";
	next;
    }
    if($opt{v}){
      ##QC step 6
      #skip redundant variants with incorrect transcript
      if($target < 0){
	unless(index($consq, "SYNONYMOUS_CODING") != -1) {
	  print FUNCE "TRANSCRIPT DOES NOT MATCH DATA; noncoding annotation inside exon \t$variant\n";
	  next;
	}	
      }
      if($target >= 0){
	if(index($consq, "SYNONYMOUS_CODING") != -1) {
	  print FUNCE "TRANSCRIPT DOES NOT MATCH DATA; coding annotation outside exon \t$variant\n";
	  next;
	}	
      }
    }
   
    ###########################################################    
    if($var eq $known){
	unless($cdna eq "." or $cdna eq ""){push(@totcdna,$cdna) unless ( $cdna ~~ @totcdna );}
	unless($prot eq "." or $prot eq ""){push(@totprot,$prot) unless ( $prot ~~ @totprot );}
	unless($cdna eq "." or $cdna eq ""){push(@totenst,$enst) unless ( $enst ~~ @totenst );}
	unless($consq eq "." or $consq eq ""){push(@totconsq,$consq) unless ( $consq ~~ @totconsq );}
    }
    else{
	$known = $var;
	@totcdna=();
	@totprot=();
	@totenst=();
	@totconsq=();
	unless($cdna eq "." or $cdna eq ""){push(@totcdna,$cdna);}
	unless($prot eq "." or $prot eq ""){push(@totprot,$prot);}
	unless($cdna eq "." or $cdna eq ""){push(@totenst,$enst);}
	unless($consq eq "." or $consq eq ""){push(@totconsq,$consq);}
    }
    my $cdnas=join(",",@totcdna);
    my $prots=join(",",@totprot);
    my $ensts=join(",",@totenst);
    my $consqs=join(",",@totconsq);

    #annotations
    $an1{$var}="HGNC=$hgnc" unless $hgnc eq "." or $hgnc eq ""; 
    $an2{$var}="CDNA=$cdnas" unless $cdna eq "." or $cdna eq "";
    $an3{$var}="VARAA=$prots" unless $prot eq "." or $prot eq "";
    $an4{$var}="TYPE=$type" unless $type eq "." or $type eq "";
    $an5{$var}="CONSQ=$consqs" unless $consq eq "." or $consq eq "";
    #$an5{$var}="ET=$ensts" unless $prot eq "." or $prot eq "";
    $an6{$var}="ET=$ensts" unless $cdna eq "." or $cdna eq "";

    #save clean func file
    print FUNCN "$variant\n";

    ##skip correct redundant lines for processing and statistics: 
    #	multiple alleles printed separately with same genotype info.
    #	for each alleles multiple copies with different ENSEMBLT 
    unless($genotype{$var}){
      ##interpret caller info for quality control, count stats
      my $GQ = my $biq = my $biqi = "NA";
      my $baseQ = my $QD = my $MQ = my $FS = my $HS = my $MQRS = my $RPRS = my $VQSLOD = my $IC="NA";
      my $NF = my $NR = my $NRS = my $NFS = my $HP = my $DGQ = "NA";
      my $VDB = my $MMQ = my $FQ = my $SP ="NA";
      my @PV=("NA","NA","NA","NA");
      my @PL=("NA","NA","NA");
      foreach my $c(@callers){
	$c =~ s/;\[//;#remove caller info brackets 
	$c =~ s/\[//;
	my (%stats);
	my @caller=split(":",$c);
	my $calln=shift(@caller);
	my $rest=join(':',@caller);
	my @callerinfo=split(";",$rest);
	foreach my $field(@callerinfo){
		my @pair=split("=",$field);
        	$stats{$pair[0]}=$pair[1];
	}
	my $locorg=$stats{"LOC_ORG"} if $stats{"LOC_ORG"};
	#go through callers	
 	if($calln eq "GATK" or ($calln eq "GATK_SANGER" and $cn==1)){
	   @gloc=split(",",$locorg);# original locus info
	   $isgatk = 1;
	   $GQ=$stats{"GQ"} if $stats{"GQ"};
	   $baseQ=$stats{BaseQRankSum} if $stats{BaseQRankSum};
	   $QD=$stats{QD} if $stats{QD};
	   $MQ=$stats{MQ} if $stats{MQ};
	   $FS=$stats{FS} if $stats{FS};
	   $HS=$stats{HaplotypeScore} if $stats{HaplotypeScore};
	   $MQRS=$stats{MQRankSum} if $stats{MQRankSum};
	   $RPRS=$stats{ReadPosRankSum} if $stats{ReadPosRankSum};
	   $VQSLOD=$stats{VQSLOD} if $stats{VQSLOD};
	   $IC=$stats{InbreedingCoeff} if $stats{InbreedingCoeff}; 
	   $filter=$stats{FILTER} if $stats{FILTER};
	   $gatk++;
	   if($type eq "SNP"){
		$biq = 0;
		if($stats{QD} and $stats{QD} > 2.0){
			$biq++;
			$qdpass++;
		}
		elsif($stats{QD} and $stats{QD} < 2.0){
			$qdfail++;
		}
		if($stats{MQ} and $stats{MQ} > 40.0){
			$biq++;
			$mqpass++;
		}
		elsif($stats{MQ} and $stats{MQ} < 40.0){
			$mqfail++;
		}
		if($stats{FS} and $stats{FS} < 60.0){
			$biq++;
			$fspass++;
		}
		elsif($stats{FS} and $stats{FS} > 60.0){
			$fsfail++;
		}
		if($stats{HaplotypeScore} and $stats{HaplotypeScore} < 13.0){
			$biq++;
			$hspass++;
		}
		if($stats{HaplotypeScore} and $stats{HaplotypeScore} > 13.0){
			$hsfail++;
		}
		if($stats{MQRankSum} and $stats{MQRankSum} > -12.5){
			$biq++;
			$mqrspass++;
		}
		if($stats{MQRankSum} and $stats{MQRankSum} < -12.5){
			$mqrsfail++;
		}
		if($stats{ReadPosRankSum} and $stats{ReadPosRankSum} > -8.0){
			$biq++;
			$rprspass++;
		}
		if($stats{ReadPosRankSum} and $stats{ReadPosRankSum} < -8.0){
			$rprsfail++;
		}
		if($stats{VQSLOD} and $stats{VQSLOD} >= 3){
			$vqslodpass++;
		}
		if($stats{VQSLOD} and $stats{VQSLOD} < 3){
			$vqslodfail++;
		}
	   }
	   elsif($type eq "DEL" or $type eq "INS" or $type eq "INDEL"){
		$biqi = 0;
		if($stats{QD} and $stats{QD} > 2.0){
			$biqi++;
			$iqdpass++;
		}
		if($stats{QD} and $stats{QD} < 2.0){
			$iqdfail++;
		}
		if($stats{FS} and $stats{FS} < 200.0){
			$biqi++;
			$ifspass++;
		}
		if($stats{FS} and $stats{FS} > 200.0){
			$ifsfail++;
		}
		if($stats{ReadPosRankSum} and $stats{ReadPosRankSum} > -20.0){
			$biqi++;
			$irprspass++;
		}
		if($stats{ReadPosRankSum} and $stats{ReadPosRankSum} < -20.0){
			$irprsfail++;
		}
		if($stats{InbreedingCoeff} and $stats{InbreedingCoeff} > -0.8){
			$biqi++;
			$icpass++;
		}
		if($stats{InbreedingCoeff} and $stats{InbreedingCoeff} < -0.8){
			$icfail++;
		}
		if($stats{VQSLOD} and $stats{VQSLOD} >= 3){
			$ivqslodpass++;
		}
		if($stats{VQSLOD} and $stats{VQSLOD} < 3){
			$ivqslodfail++;
		}
	   }
	   if($cn==1){
		$gsolo++;
		printf GATK2 "$funcf\t$calln\t$type\t$cover\t$refF\t$altF\t$GQ\t$QD\t$MQ\t$FS\t$HS\t$MQRS\t$RPRS\t$baseQ\t$VQSLOD\t$IC\n";
	   }
	   else{printf GATK "$funcf\t$calln\t$type\t$cover\t$refF\t$altF\t$GQ\t$QD\t$MQ\t$FS\t$HS\t$MQRS\t$RPRS\t$baseQ\t$VQSLOD\t$IC\n";}
	}
	elsif($calln eq "MPILEUP"){
	   @mloc=split(",",$locorg);# original locus info
	   $ismpileup = 1;
	   $VDB=$stats{VDB} if $stats{VDB};
	   $MMQ=$stats{MQ} if $stats{MQ};
	   $FQ=$stats{FQ} if $stats{FQ};
	   if($stats{PV4}){
		my $PV4=$stats{PV4};
		@PV=split(/,/,$PV4);
	   }
	   if($stats{PL}){
		my $PL3=$stats{PL};
		@PL=split(/,/,$PL3);
	   }
	   $SP=$stats{SP} if $stats{SP};
	   $mpileup++;
	   if($cn==1){
		$msolo++;
		printf MPILEUP2 "$funcf\t$calln\t$type\t$cover\t$refF\t$altF\t$MMQ\t$FQ\t$PV[0]\t$PV[1]\t$PV[2]\t$PV[3]\t$PL[0]\t$PL[1]\t$PL[2]\t$SP\n";
	   }
	   else{printf MPILEUP "$funcf\t$calln\t$type\t$cover\t$refF\t$altF\t$MMQ\t$FQ\t$PV[0]\t$PV[1]\t$PV[2]\t$PV[3]\t$PL[0]\t$PL[1]\t$PL[2]\t$SP\n";}
	}
	elsif($calln eq "DINDEL" or $calln eq "DINDEL_SANGER"){
	   @dloc=split(",",$locorg);# original locus info
	   $isdindel = 1;
	   $NF=$stats{NF} if $stats{NF};
	   $NR=$stats{NR} if $stats{NR};
	   $NRS=$stats{NRS} if $stats{NRS};
	   $NFS=$stats{NFS} if $stats{NFS};
	   $HP=$stats{HP} if $stats{HP};
	   $DGQ=$stats{GQ} if $stats{GQ};
	   $dindel++;
	   if($cn==1){
		$dsolo++; 
		printf DINDEL2 "$funcf\t$calln\t$type\t$cover\t$refF\t$altF\t$DGQ\t$NF\t$NR\t$NRS\t$NFS\t$HP\n";
	   }		
	   else{printf DINDEL "$funcf\t$calln\t$type\t$cover\t$refF\t$altF\t$DGQ\t$NF\t$NR\t$NRS\t$NFS\t$HP\n";}
	}
	elsif($calln eq "MRSSVSP"){
	   $ismrssvsp = 1;
	   $mrssvsp++;
	   if($cn==1){$mrsolo++;}
	}
      }
      printf ALL "$funcf\t$type\t$isgatk\t$ismpileup\t$isdindel\t$cover\t$refF\t$altF\t$GQ\t$QD\t$MQ\t$FS\t$HS\t$MQRS\t$RPRS\t$baseQ\t$VQSLOD\t$IC\t$MMQ\t$FQ\t$PV[0]\t$PV[1]\t$PV[2]\t$PV[3]\t$PL[0]\t$PL[1]\t$PL[2]\t$SP\t$DGQ\t$NF\t$NR\t$NRS\t$NFS\t$HP\n";

      #count filter passes 
      if($biq eq 0){$gatk0++;}
      elsif($biq eq 1){$gatk1++;}
      elsif($biq eq 2){$gatk2++;}
      elsif($biq eq 3){$gatk3++;}
      elsif($biq eq 4){$gatk4++;}
      elsif($biq eq 5){$gatk5++;}
      elsif($biq eq 6){$gatk6++;}
      if($biqi eq 0){$igatk0++;}
      elsif($biqi eq 1){$igatk1++;}
      elsif($biqi eq 2){$igatk2++;}
      elsif($biqi eq 3){$igatk3++;}
      elsif($biqi eq 4){$igatk4++;}

      ##define hashes that stores all fixed data
      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER
      $fixed1{$var}=$chrom."\t".$opos."\t".$id."\t".$oref."\t".$oalt."\t".$qual."\t".$filter unless $fixed1{$var};

      ##GENOTYPE
      #get allele counts	
      my ($sAC1,$sAC2);
      if($data[13] eq "1|1" or $data[13] eq "1/1"){
	if($is_male == 1 and ($chrom eq "X" or $chrom eq "Y")){$sAC1=1;}
	else{$sAC1=2;}
      }
      elsif($data[13] eq "1|2" or $data[13] eq "1/2" or $data[13] eq "2|1"){
	$sAC1=1;
	$sAC2=1;
      }
      elsif($data[13] eq "0|1" or $data[13] eq "0/1" or $data[13] eq "1|0"){$sAC1=1;}
      elsif($data[13] eq "0|2" or $data[13] eq "0/2" or $data[13] eq "2|0"){
	$sAC1=0;
	$sAC2=1;
      }
      elsif($data[13] eq "2|2" or  $data[13] eq "2/2"){
	if($is_male == 1 and ($chrom eq "X" or $chrom eq "Y")){
	  $sAC1=0;
	  $sAC2=1;
	}
	else{
	  $sAC1=0;
	  $sAC2=2;
	}
      }
      else{print FUNCE "OTHER_ERROR\t$variant\n";}
      #define format field
      $format="GT:GQ:VQ:CP:BF:NC:DP";
      #count only non-redundant variants
      #unless($genotype{$var}){
      #allele counts
      if(defined $sAC1){#total number of first allele counts
	if($AC1{$var}){$AC1{$var}=$AC1{$var}+$sAC1;}
	else{$AC1{$var}=$sAC1;}
      }
      if(defined $sAC2){#total number of second allele counts
	if($AC2{$var}){$AC2{$var}=$AC2{$var}+$sAC2;}
	else{$AC2{$var}=$sAC2;}
      }
      $alts{$var}=$nalt;
      #genotype
      $genotype{$var}=$data[13];
      my $cp = 0;
      if($biq eq "NA"){$cp = $biqi;}
      else{$cp = $biq;}
      #genotype field includes: GQ=genotyping quality, number of BRoad Institute criteria fulfilled, frequency of reference allele, frequency of alternative allele and coverage
      my $gtpfield=$genotype{$var}.":".$GQ.":".$VQSLOD.":".$cp.":".nearest(.01, $altF).":".$cn.":".$cover;
      $genotyping{$var}[$NS] = $gtpfield;

    }#end of processing of unique variant (first annotation corresponding to protein coding ENST)

   }
  }#end of data set (func file) loop

  ###########################################################
  ##variant stats
  if($is_male == 1){
    $males++;
    $sex[$NS]="M";
  }
  else{
    $sex[$NS]="F";
  }
  close (FUNC);
  close (FUNCE);
  close (FUNCN);
  my $total=$gatk0+$gatk1+$gatk2+$gatk3+$gatk4+$gatk5+$gatk6;
  printf GSTAT "$funcf\t$total\t$gatk\t$mpileup\t$dindel\t$mrssvsp\t$gsolo\t$msolo\t$dsolo\t$mrsolo\n";
  printf GSNP "$funcf\t$gatk0\t$gatk1\t$gatk2\t$gatk3\t$gatk4\t$gatk5\t$gatk6\t$qdpass\t$qdfail\t$mqpass\t$mqfail\t$fspass\t$fsfail\t$hspass\t$hsfail\t$mqrspass\t$mqrsfail\t$rprspass\t$rprsfail\t$vqslodpass\t$vqslodfail\n";
  printf GIND "$funcf\t$igatk0\t$igatk1\t$igatk2\t$igatk3\t$igatk4\t$iqdpass\t$iqdfail\t$ifspass\t$ifsfail\t$irprspass\t$irprsfail\t$ivqslodpass\t$ivqslodfail\t$icpass\t$icfail\n";
}
close (GATK);
close (GATK2);
close (MPILEUP);
close (MPILEUP2);
close (DINDEL);
close (DINDEL2);
close (GSTAT);
close (GSNP);
close (GIND);


###########################################################
##sumarize all
print GOUT "\n";#finish header
#summarize all data sets
my $AN0=$NS*2;
my $ANy=$males;
my $ANx=($NS-$males)*2+$males;
#
for my $var (sort keys %fixed1){
  my ($info,$AN);
  my @temp=split("-",$var);
  my $chr=$temp[0];
  if($chr eq "Y"){$AN=$ANy;}
  elsif($chr eq "X"){$AN=$ANx;}
  else{$AN=$AN0;}
  #info: NS=Number of samples with data; AN=total number of alleles in called genotypes; AC=allele count in genotypes, for each ALT allele, in the same order as listed; AF=allele frequency for each ALT allele in the same order as listed
  if($AC2{$var}){
    my $AF1=nearest(.0001, $AC1{$var}/$AN);
    my $AF2=nearest(.0001, $AC2{$var}/$AN);
    my $AF="$AF1,$AF2";
    my $AC="$AC1{$var},$AC2{$var}";
    if($alts{$var} == 3){$AF=$AF.",".0;}
    if($alts{$var} == 3){$AC=$AC.",".0;}
    $info="NS=$NS;AN=$AN;AC=$AC;AF=$AF";
    if($an1{$var}){$info=$info.";".$an1{$var};}
    if($an2{$var}){$info=$info.";".$an2{$var};}
    if($an3{$var}){$info=$info.";".$an3{$var};}
    if($an4{$var}){$info=$info.";".$an4{$var};}
    if($an5{$var}){$info=$info.";".$an5{$var};}
    if($an6{$var}){$info=$info.";".$an6{$var};}
  }
  else{
    my $AF1=nearest(.0001, $AC1{$var}/$AN);
    $info="NS=$NS;AN=$AN;AC=$AC1{$var};AF=$AF1";
    if($an1{$var}){$info=$info.";".$an1{$var};}
    if($an2{$var}){$info=$info.";".$an2{$var};}
    if($an3{$var}){$info=$info.";".$an3{$var};}
    if($an4{$var}){$info=$info.";".$an4{$var};}
    if($an5{$var}){$info=$info.";".$an5{$var};}
    if($an6{$var}){$info=$info.";".$an6{$var};}
  }
  printf GOUT  "$fixed1{$var}\t$info\t$format";
  for(my $i=1;$i<=$NS;$i++){
    unless ($genotyping{$var}[$i]){
      if($sex[$NS] eq "M"){#correct genotype of haploid alleles
	if($chr eq "X" or $chr eq "Y"){$genotyping{$var}[$i]="0";}
	else{$genotyping{$var}[$i]="0/0:.:.:.:.:.:.";}
      }
      else{
	$genotyping{$var}[$i]="0/0:.:.:.:.:.:.";
      }
    }
    else{
      if($sex[$NS] eq "M"){#correct genotype of haploid alleles
	if($chr eq "X" or $chr eq "Y"){
	  my @gt=split(":",$genotyping{$var}[$i]);
	  if($gt[0] eq "1|1" or $gt[0] eq "1/1" or $gt[0] eq "0|1" or $gt[0] eq "0/1" or $gt[0] eq "1|0"){$gt[0]=1;}
	  if($gt[0] eq "2|2" or $gt[0] eq "2/2" or $gt[0] eq "0|2" or $gt[0] eq "0/2" or $gt[0] eq "2|0"){$gt[0]=2;}
	  $genotyping{$var}[$i]=join(":",@gt);
	}
      }   
    }
    printf GOUT "\t$genotyping{$var}[$i]";
  }
  print GOUT "\n";
}

sub usage(){
print STDERR << "EOF";
Usage: func2vcf.pl -f [file with list of func files] -o [name of output file]
 -h     : help message
 -t	: remove variant outside of targets by more than 100bp
 -v	: takes only transcription variants matching the targeted region

EOF
exit;
}
sub tutej_any { ( grep $_, @_ ) < 2 }

