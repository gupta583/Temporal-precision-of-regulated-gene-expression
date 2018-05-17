#!/bin/bash

#rm seed.dat
for i in {1..50}
do
gfortran reaction.f90
echo -733$RANDOM | ./a.out > data/file_$i.dat

done

