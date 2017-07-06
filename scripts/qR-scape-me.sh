#!/bin/bash

wraprscape=$HOME/Perl/MyModules/Bio-RNA-RNAaliSplit/scripts/wrap_R-scape.pl
outdir=/scr/hell/mtw/Work/Nitrososphaera/aln
stat=RAFS
log=/scr/hell/mtw/Work/Nitrososphaera/aln/rscape.log

for file in $(find . -name \*stk -print0 | xargs -0) 
do	
	 qsub perl $wraprscape -a $file -o $outdir -s $stat -l $log --nofigures
done

