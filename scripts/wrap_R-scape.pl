#!/usr/bin/env perl
# Last changed Time-stamp: <2017-07-06 19:52:14 mtw>
# -*-CPerl-*-
#
# A wrapper for R-scape
#
# usage: wrap_R-scape.pl -a myaln.stk --statistics RAFS
#

use version; our $VERSION = qv('0.05.2');
use strict;
use warnings;
use Bio::RNA::RNAaliSplit::WrapRscape;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Path::Class;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $show_version = 0;
my ($id,$stkfile,$logfile);
my $handle=undef;
my $nofigures = 0;
my $outdir = "rscape";
my $stat = "GTp";
my $have_logfile=0;
use diagnostics;

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("a|aln=s"        => \$stkfile,
					   "s|statistic=s"  => \$stat,
					   "o|out=s"        => \$outdir,
					   "l|log=s"        => \&set_logfile,
					   "nofigures"      => sub{$nofigures = 1},
                                           "version"        => sub{$show_version = 1},
                                           "man"            => sub{pod2usage(-verbose => 2)},
                                           "help|h"         => sub{pod2usage(1)}
                                           );

if ($show_version == 1){
  print "wrap_R-scape.pl $VERSION\n";
  exit(0);
}

sub set_logfile {
  $logfile=$_[1];
  $have_logfile=1;
}

if ($have_logfile == 1){
  open($handle, ">>", $logfile)
    || die "$0: can't open $logfile for appending: $!";
}
else{
  open($handle, ">&", STDOUT)
    || die "$0: can't open handle for STDOUT: $!";
}


my @odir = split /\//, $outdir;

unless (-f $stkfile){
  warn "Could not find input file provided via -a|-aln option";
  pod2usage(-verbose => 0);
}

my $r = Bio::RNA::RNAaliSplit::WrapRscape->new(ifile => $stkfile,
					       statistic => $stat,
					       nofigures => $nofigures,
					       odir => \@odir,
					      );
#print Dumper($r);

if ($r->cseq <= 1){ # stk file had only one sequence
  print $handle join "\t", ($stkfile, "n/a");
}
else{ # normal stk file
  my $ratio;
  $r->nbpairs > 0 ? ($ratio = $r->TP/$r->nbpairs) : ($ratio = 0);
  print $handle join "\t", ($stkfile,$r->statistic,$r->TP,$r->alen,$r->nbpairs,$r->nseq, $ratio);
}
print $handle "\n";

__END__

=head1 NAME

wrap_R-scape.pl - Wrapper script for R-scape >v0.6.1

=head1 SYNOPSIS

wrap_R-scape.pl [--aln|-a I<FILE>] [--statistic|-s I<STRING>] [--nofigures] [options]

=head1 DESCRIPTION

This is a Perl wrapper for R-scape >v0.6.1
(http://eddylab.org/R-scape/). It accepts a single multiple sequence
alignment in Stockholm format and runs R-scape on it, determining
statistically significant covarying base pairs (SSCBP).

The applied covariance statistic, the number of SSCBP, the lenth of
the alignment, as well as the number of base pairs and sequences is
reopted on STDOUT. In addition the fraction of SSCBP over all base
pairs and the input file name is listed for further processing. All
R-scape output files are written to a user-defined directory.

This script was inteded as simple R-scape warpper. As such, it does
not implement (and pass through) all R-scape options.

Please ensure that R-scape >v0.6.1 is installed on your machine and
available for your Perl interpreter.

=head1 OPTIONS

=over

=item B<--aln|-a>

A multiple sequence alignment in Stockholm format.

=item B<--statistic|-s>

The covariation statistic used by R-scape. Allowed values are: 'GT',
'MI', 'MIr', 'MIg', 'CHI', 'OMES', 'RAF', 'RAFS'. Appending either 'p'
or 'a' to any of them calculates its average product correction and
average sum correction, respctively (e.g. GTp or GTa). See the R-scape
manual for details.

=item B<--nofigures>

Turn off production of graphical R-scape output

=item B<--log|-l>

Specify log file; output is written to STDOUT unless specified

=item B<--version>

Show wrap-R-scape.pl version and exit

=back

=head1 SEE ALSO

The L<R-scape
manual|http://eddylab.org/software/rscape/R-scape_userguide.pdf> and
the L<R-scape
paper|http://eddylab.org/publications/RivasEddy16/RivasEddy16-preprint.pdf>.

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt> and E<lt>michael.wolfinger@univie.ac.atE<gt>

=cut
