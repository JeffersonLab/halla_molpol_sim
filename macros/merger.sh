#!/bin/bash
#
dir=$HOME/gitdir/qsim/build/output/
i=2
while [ $i -le 50 ];
do
	hadd $dir/result.root $dir/out$i/qsim_output.root $dir/Result.root
	let i=i+1
	mv $dir/result.root $dir/Result.root
done
