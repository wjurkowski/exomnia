#!/usr/bin/perl -w
use strict;
use warnings;

#if ($#ARGV != 2) {die "Program used with parameters [pdb] [list of residues] [SIFTS mapping file]\n";}
if ($#ARGV != 3) {die "Program used with parameters [pdb] [chid] [list of residues] [SIFTS mapping file]\n";}

#uniprotAC
my $pdb=$ARGV[0];
my $chid=$ARGV[1];

#read list of residues 
my $vars=$ARGV[2]; 
my @snv=open_file($vars);

#parse PDB swissport mapping
my @sifts=open_file($ARGV[3]);
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
print "cycy $mypdbch $pdb2uni{$mypdbch}\n";
my @mapping=split(/\t/,$sifts[$pdb2uni{$mypdbch}]);
my $UPseq_s=$mapping[7];
my $PDBseq_s=$mapping[5];
my $shift=$PDBseq_s-$UPseq_s;

my ($string);
foreach my $lin(@snv){
	my @aa=split(/\t/,$lin);
	my $pos=$aa[1]+$shift;
	if ($string){
	  $string=$string."+".$pos unless $pos < 1;
	}
	else{
	  $string=$pos unless $pos < 1;
	}
}

#prepare pymol script
open(PML, ">mapListofResidues.pml") or die "Can not open an output file: $!";
print PML "\n";
print PML "log_open mapListofResidues.log\n";
print PML "load pdb$pdb.ent\n";
print PML "show ribbon\n";
print PML "hide lines\n";
#	color groups of aa
#	group 1
print PML "set_color dblue, [30,144,255]\n";#dodger blue
print PML "set_color fgreen, [34,139,34]\n";#forest green
print PML "set_color fire, [178,34,34]\n";#firebrick
print PML "set_color orange, [255,165,0]\n";#orange
print PML "set_color orchid, [153,50,204]\n";#dark orchid
print PML "set_color vred, [208,32,144]\n";#violet red
print PML "select $ARGV[2],resid $string\n";
print PML "color dblue, $ARGV[2]\n";
print PML "show sticks, vars\n";

sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}



