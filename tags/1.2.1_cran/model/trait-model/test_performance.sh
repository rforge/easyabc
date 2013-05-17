#!/bin/bash

echo "test of model performance" >time_performance.txt
date >> time_performance.txt


size=3


# I/O one by one
for ((i=0 ; i<$size ; i++))
do
nbsimu=$((10**$i))
echo $nbsimu

begintime=`date +%s`
for ((j=0 ; j<$nbsimu ; j++))
do
./trait_model
done
endtime=`date +%s`
echo "ioalone;$nbsimu;$begintime;$endtime;$(($endtime-$begintime))">> time_performance.txt


done




# I/O by group
for ((i=0 ; i<$size ; i++))
do
nbsimu=$((10**$i))
echo $nbsimu

echo "$nbsimu">input2
for ((j=0 ; j<$nbsimu ; j++))
do
cat input2_0>>input2
done
cp "input_ref_$nbsimu" input_ref
begintime=`date +%s`
./trait_model2
endtime=`date +%s`
echo "iogroup;$nbsimu;$begintime;$endtime;$(($endtime-$begintime))">> time_performance.txt


done





# reference
for ((i=0 ; i<$size ; i++))
do
nbsimu=$((10**$i))
echo $nbsimu

cp "input_ref_$nbsimu" input_ref
begintime=`date +%s`
./trait_model_ref
endtime=`date +%s`
echo "ref;$nbsimu;$begintime;$endtime;$(($endtime-$begintime))">> time_performance.txt


done


# todo test from R


