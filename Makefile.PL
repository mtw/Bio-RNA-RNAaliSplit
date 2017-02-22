use 5.010;
use lib 'inc';
use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share dist => 'share';

my @scripts = qw (RNAalisplit.pl);

WriteMakefile(
	      NAME               => 'Bio::RNA::AliSplit',
	      AUTHOR             => q{Michael T. Wolfinger <michael@wolfinger.eu>},
	      VERSION_FROM       => 'lib/Bio/RNA/AliSplit.pm',
	      ABSTRACT_FROM      => 'lib/Bio/RNA/AliSplit.pm',
	      LICENSE            => 'perl_5',
	      PL_FILES           => {},
	      MIN_PERL_VERSION   => 5.010,
	      CONFIGURE_REQUIRES => {
				     'ExtUtils::MakeMaker' => 0,
				    },
	      BUILD_REQUIRES     => {
				     'Test::More' => 0,
				    },
	      PREREQ_PM          => {
				     "Getopt::Long"	     => '2.0',
				     "Bio::Root::Version"    => '1.00690001',
				    },
	      'EXE_FILES'        => [ map { "scripts/$_" } @scripts ],
	      dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
	      clean => { FILES => 'Bio-RNA-AliSplit-*' },
	     );

package MY;
use File::ShareDir::Install qw(postamble);
