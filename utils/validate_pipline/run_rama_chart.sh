#!/bin/bash

#Script to run the missing phenix.striage validations
#runs if no xtriage log exist and structure has an mtz.
#Please, run this from the utils folder
#x is the full path to the folder
#c is the pdb code
#e is the entry in lower case (for pdb_redo)

cd ../../pdb

base=$(pwd)
fivedaysago=$(date -d 'now - 5 days' +%s)

for x in */*/* ;
do
    c=`basename $x`;
    echo $x
    #e="${c,,}" ;
    cd $x ;
    if [[ ( "$c" == "7mhl" ) || ( "$c" == "7mhm" ) || ( "$c" == "7mhq" ) || ( "$c" == "6zp5" ) || 
    	( "$x" == *"txt"* ) || ( "$x" ==  *"xlsx"* ) || ( "$x" == *"fasta"* ) ]]; then
    	cd $base
    	continue
    fi
    file_time=$(date -r "validation/molprobity/$c.H.pdb" +%s)
    if [[ "$file_time" > "$fivedaysago" ]] ; then
  		#echo "$c has been processed"
  		cd $base
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