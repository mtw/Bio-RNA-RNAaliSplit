#!/usr/bin/env perl
# Last changed Time-stamp: <2017-03-21 04:40:52 michl>
# -*-CPerl-*-
#
# usage: split_stockholm.pl -a myln.stk
#

use version; our $VERSION = qv('0.05_01');
use strict;
use warnings;
use Bio::RNA::RNAaliSplit::WrapRNAalifold;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Path::Class;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $format = "ClustalW";
my $show_version = 0;
my ($id,$alifile);
my $nr = 1;
my %stks=();

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("aln|a=s"    => \$alifile,
					   "version"    => sub{$show_version = 1},
                                           "man"        => sub{pod2usage(-verbose => 2)},
                                           "help|h"     => sub{pod2usage(1)}
                                           );

if ($show_version == 1){
  print "split_stockholm $VERSION\n";
  exit(0);
}

unless (-f $alifile){
  warn "Could not find input file provided via --aln|-a option";
  pod2usage(-verbose => 0);
}

open my $input, "<", $alifile or die $!;
my $fn = "$nr.aln";

# parse multi-Stockholm file into %stks hash
while(<$input>){
  chomp;
  my $line = $_;
  unless ($line eq "//"){
    if ($line =~ /\#=GF\sID\s(.+)/){
      $id = $1;
 #     print "+++++++++++ $id +++++++++++\n";
    }
    push @{$stks{$nr}{data}}, $line;
    $stks{$nr}{id}=$id;
    #print "$line\n";
  }
  else{
    $nr++;
 #   print "---------------------------------\n";
  }
}
close $input;

# process %stks hash and write individual Stockholm files
foreach my $i (keys %stks){
  my $name = $stks{$i}{id}.".stk";
  open my $alnfile, ">", $name or die $!;
  foreach my $l ( @{$stks{$i}->{data}} ){
    print $alnfile $l."\n";
  }
  print $alnfile "//\n";
  #print Dumper (\@{$stks{$i}->{data}});
  close $alnfile;
  #print "----\n";
}


