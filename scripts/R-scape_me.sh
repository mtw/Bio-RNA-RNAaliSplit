#!/bin/bash

wraprscape=$HOME/Perl/MyModules/Bio-RNA-RNAaliSplit/scripts/wrap_R-scape.pl
stks=$(find . -name *stk)

for file in $stks
do	
	$wraprscape -a $file -o test -s RAFS --nofigures
done

