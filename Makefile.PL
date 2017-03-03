use 5.010;
use lib 'inc';
use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share 'share';

my @scripts = qw (RNAalisplit.pl);

WriteMakefile(
	      NAME               => 'Bio::RNA::RNAaliSplit',
	      AUTHOR             => q{Michael T. Wolfinger <michael@wolfinger.eu>},
	      VERSION_FROM       => 'lib/Bio/RNA/RNAaliSplit.pm',
	      ABSTRACT_FROM      => 'lib/Bio/RNA/RNAaliSplit.pm',
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
				     "Array::Set"=> '0',
				     "Bio::Root::Version" => '1.00690001',
				     "Bio::AlignIO" => '0',
				     "Carp" => '0',
				     "Data::Dumper" => '0',,
				     "Digest::MD5" => '0',
				     "File::Basename" => '0',
				     "File::Path" => '0',
				     "File::Spec" => '0',
				     "Getopt::Long" => '2.0',
				     "IO::Dir" => '0',
				     "IPC::Cmd" => '0',
				     "Moose" => '0',
				     "Moose::Role" => '0',
				     "Moose::Util::TypeConstraints" => '0',
				     "Path::Class" => '0',
				     "Path::Class::Dir" => '0',
				     "Path::Class::File" => '0',
				     "RNA" => "2.3.3",
				     "Storable" => '0',
				    },
	      'EXE_FILES'        => [ map { "scripts/$_" } @scripts ],
	      dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
	      clean => { FILES => 'Bio-RNA-AliSplit-*' },
	     );

package MY;
use File::ShareDir::Install qw(postamble);