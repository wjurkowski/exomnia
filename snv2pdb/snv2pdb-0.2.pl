#!/usr/bin/perl -w
use strict;
use warnings;

#if ($#ARGV != 2) {die "Program used with parameters [pdb] [list of residues] [SIFTS mapping file]\n";}
if ($#ARGV != 4) {die "Program used with parameters [pdb] [chid] [rare SNV] [common SNV] [SIFTS mapping file]\n";}

#uniprotAC
my $pdb=$ARGV[0];
my $chid=$ARGV[1];

#read list of rare variants
my @rare=open_file($ARGV[2]);
#read list of common variants
my @common=open_file($ARGV[3]);

#parse PDB swissport mapping
my @sifts=open_file($ARGV[4]);
my (%pdb2uni,%uniAC2pdb);
for (my $i=0;$i<=$#sifts;$i++){
	my @lin=split(/\t/,$sifts[$i]);
	my $pdbch=$lin[0].$lin[1];
	$pdb2uni{$pdbch}=$i;
#print "kaka $pdbch $pdb2uni{$pdbch}\n"
#	$uniAC2pdb{$lin[2]}= $pdbch;	
}

#find pdbID for query uniAC
#my $mypdbch = $uniAC2pdb{$uniAC};
my $mypdbch = $pdb.$chid;
#get mapping 
#print "cycy $mypdbch $pdb2uni{$mypdbch}\n";
my @mapping=split(/\t/,$sifts[$pdb2uni{$mypdbch}]);
my $UPseq_s=$mapping[7];
my $PDBseq_s=$mapping[5];
#print "$UPseq_s $PDBseq_s\n";
my $shift=$PDBseq_s-$UPseq_s;
print "PDB seq start - Uniprot seq start = $shift\n";

my ($resr,$resc);
foreach my $lin(@rare){
	my @aa=split(/\t/,$lin);
	my $pos=$aa[1]+$shift;
	if ($resr){
	  $resr=$resr."+".$pos unless $pos < 1;
	}
	else{
	  $resr=$pos unless $pos < 1;
	}
}
foreach my $lin(@common){
	my @aa=split(/\t/,$lin);
	my $pos=$aa[1]+$shift;
	if ($resc){
	  $resc=$resc."+".$pos unless $pos < 1;
	}
	else{
	  $resc=$pos unless $pos < 1;
	}
}
unless ($resc or $resr){print "none of SNVs present in pdb structure\n";}
#prepare pymol script
open(PML, ">displaySNVs.pml") or die "Can not open an output file: $!";
#print PML "log_open mapListofResidues.log\n";
#print PML "load pdb$pdb.ent\n";
#	define new colors
print PML "set_color dblue, [30,144,255]\n";#dodger blue
print PML "set_color fgreen, [34,139,34]\n";#forest green
print PML "set_color fire, [178,34,34]\n";#firebrick
print PML "set_color orange, [255,165,0]\n";#orange
print PML "set_color orchid, [153,50,204]\n";#dark orchid
print PML "set_color vred, [208,32,144]\n";#violet red
#fetch biological assembly
print PML "fetch $pdb, type=pdb1, async=0\n";
print PML "split_state $pdb\n";
#modify general display
print PML "show cartoon\n";
print PML "hide lines\n";
print PML "remove solvent\n";
print PML "util.cbc\n";
print PML "show spheres, organic\n";
print PML "util.cbag organic\n";
#display rare and common variants
print PML "select rare,resid $resr\n" if ($resr);
print PML "select common,resid $resc\n" if ($resc);;
#	group 1
print PML "color orchid, rare\n";
print PML "show sticks, rare\n";
#	group 2
print PML "color orange, common\n";
print PML "show sticks, common\n";


sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}



