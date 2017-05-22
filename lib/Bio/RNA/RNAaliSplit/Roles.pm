# -*-CPerl-*-
# Last changed Time-stamp: <2017-05-22 23:43:51 mtw>

package Bio::RNA::RNAaliSplit::Roles;

use version; our $VERSION = qv('0.05_02');
use Moose::Util::TypeConstraints;
use Moose::Role;
use Path::Class::Dir;
use namespace::autoclean;

has 'dirnam' => ( # custom output dir name
		  is => 'rw',
		  isa => 'Path::Class::Dir',
		  predicate => 'has_dirname',
		 );

1;
