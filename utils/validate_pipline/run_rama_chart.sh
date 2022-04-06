#!/bin/bash

#Script to run the missing phenix.striage validations
#runs if no xtriage log exist and structure has an mtz.
#Please, run this from the utils folder
#x is the full path to the folder
#c is the pdb code
#e is the entry in lower case (for pdb_redo)

cd ../../pdb

base=$(pwd)


for x in */*/* ;
do
    c=`basename $x`;
    echo $x
    #e="${c,,}" ;
    cd $x ;
    if [ "$c" == "7mhl" ] || [ "$c" == "7mhm" ] || [ "$c" == "7mhq" ] || [ "$c" == "6zp5" ] ; then
    	continue
    fi
    if [[ $(find "$c" -mtime +2 -print) ]]; then
  		echo "$c has been processed"
  		continue
	fi
#	    if [ -f $c.mtz ]; then
    phenix.reduce -NOFLIP $c.pdb > validation/molprobity/$c.H.pdb
	cd validation/molprobity/
	rama_chart_pdf $c.H.pdb
	    #multichart $c.H.pdb
#            fi

cd $base
done

cd ..