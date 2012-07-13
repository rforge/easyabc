#!/bin/bash

fileName='time_performance_model.txt'
echo "test of model parthy performance" >$fileName
date >> $fileName

#../../EasyABC-R/EasyABC.R

echo "Rscript a vide"
echo "a=1">test.R
begintime=`date +%s`
Rscript test.R
endtime=`date +%s`
echo "Rvoid;0;$begintime;$endtime;$(($endtime-$begintime))">> $fileName


size=4


# I/O one by one
for ((i=0 ; i<$size ; i++))
do
nbsimu=$((10**$i))
echo $nbsimu

begintime=`date +%s`
for ((j=0 ; j<$nbsimu ; j++))
do
./parthy
done
endtime=`date +%s`
echo "ioalone;$nbsimu;$begintime;$endtime;$(($endtime-$begintime))">> $fileName


done

