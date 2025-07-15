# -*-CPerl-*-
# Last changed Time-stamp: <2025-07-15 15:39:49 mtw>

package Bio::RNA::RNAaliSplit::Subtypes;

use Moose::Util::TypeConstraints;
use Bio::AlignIO;
use Params::Coerce;
use version; our $VERSION = qv('0.12');

subtype 'Bio::RNA::RNAaliSplit::AliIO' => as class_type('Bio::AlignIO');

coerce 'Bio::RNA::RNAaliSplit::AliIO'
    => from 'HashRef'
    => via { Bio::AlignIO->new( %{ $_ } ) };
