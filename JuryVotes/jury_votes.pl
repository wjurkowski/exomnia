#!/usr/bin/perl -w
use strict;
use warnings;
if ($#ARGV != 0) {die "Program used with parameters [input file]\n";}

my @inputs=open_file($ARGV[0]);
#output file
open (GOUT, "> $ARGV[0].jury") or die "Can not create general output file: $!";


foreach my $lin (@inputs){
  my $votes=0;
  my $verdict="benign";
  my @predi=split(/\t/,$lin);
  if ($predi[6] eq "Disease") {$votes++;}
  if (($predi[4] eq "high") or ($predi[4] eq "medium")) {$votes++;}
  if (($predi[8] eq "probably damaging") or ($predi[8] eq "possibly damaging")) {$votes++;}
  if ($votes >= 2){$verdict="damaging";}
  push(@predi,$verdict);
  foreach my $val (@predi){
	printf GOUT "$val\t";
  }
  printf GOUT "\n";
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

