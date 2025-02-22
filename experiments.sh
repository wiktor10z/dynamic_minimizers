#!/bin/bash
echo $1
input_file=$1
z=$2


cd tree_minimizers
make
cd ../window_minimizers
make
cd ..

for ((l=8; l<=2048; l*=2))
do
	echo "tree" $z $l;
	./tree_minimizers/tree_minimizers -t $input_file -z $z -l $l;
done

for ((l=16; l<=2048; l*=2))
do
	echo "window" $z $l;
	./window_minimizers/window_minimizers -t $input_file -z $z -l $l;
done


cd tree_minimizers
make clean
cd ../window_minimizers
make clean
cd ..
