# -*-CPerl-*-
# Last changed Time-stamp: <2017-05-16 14:24:34 mtw>

# Bio::RNA::RNAaliSplit::FileDir.pm: A Moose Role for basic file IO

package Bio::RNA::RNAaliSplit::FileDir;

use version; our $VERSION = qv('0.05_01');
use Moose::Util::TypeConstraints;
use Moose::Role;
use Path::Class::File;
use Path::Class::Dir;

subtype 'MyFileX' => as class_type('Path::Class::File');

coerce 'MyFileX'
  => from 'Str'
  => via { Path::Class::File->new($_) };

subtype 'MyDirX' => as class_type('Path::Class::Dir');

coerce 'MyDirX'
  => from 'ArrayRef'
  => via { Path::Class::Dir->new( @{ $_ } ) };

has 'ifile' => (
		is => 'ro',
		isa => 'MyFile',
		predicate => 'has_ifile',
		coerce => 1,
	      );

has 'odir' => (
	       is => 'rw',
	       isa => 'MyDir',
	       predicate => 'has_odir',
	       coerce => 1,
	      );

has 'odirn' => ( # custom output dir name
		is => 'rw',
		isa => 'Path::Class::Dir',
		predicate => 'has_odirn',
	       );


no Moose;

1;
