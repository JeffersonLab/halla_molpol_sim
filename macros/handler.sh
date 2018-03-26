#!/bin/bash
#
location=$HOME/gitdir/qsim/build/output/
cd $location

for i in `seq 1 50`;
do
	name="out$i"
	mkdir $name
	cd $name
	qsub ../../macros/runscript.sh
	cd ..
done
