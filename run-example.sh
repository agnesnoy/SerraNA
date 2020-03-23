#!/bin/bash

#RUN extraction process.

#Remember that you should've already had compiled SerraNA (by typing make)
#and generated the elastic and structural parameters ( ./SerraNA < s_NA.in )

#Extract overalls for elastic parameters (assuming ex_NA.in hasn't been modified).
#=========================================================================================
./Extract < ex_NA.in 

#Extract elastic parameters for 4 disctinct [a,b] regions: [1,6], [7,12], [13,18], [19,24] 
#=========================================================================================
sed -i "29s/0,0/1,6/g" ex_NA.in
./Extract < ex_NA.in

sed -i "29s/1,6/7,12/g" ex_NA.in
./Extract < ex_NA.in

sed -i "29s/7,12/13,18/g" ex_NA.in
./Extract < ex_NA.in

sed -i "29s/13,18/19,24/g" ex_NA.in
./Extract < ex_NA.in

#And extract structural parameters for 2mer, 12mer and 22mer.
#=========================================================================================
sed -i "8s/   elastic_/   structural_/g" ex_NA.in
sed -i "15s/1/0/g" ex_NA.in
sed -i "29s/19,24/1,0/g" ex_NA.in
./Extract < ex_NA.in
sed -i "29s/1,0/11,0/g" ex_NA.in
./Extract < ex_NA.in
sed -i "29s/11,0/21,0/g" ex_NA.in
./Extract < ex_NA.in

#Run python script that plots 1) Bending angle for 2mer, 12mer and 22mer.
#                             2) Twist elastic constant for the 4 regions as a function of
#                                length.
#                             3) Stretch modulus as a function of length for the whole 
#                                fragment.
#
#This generates a pdf (my_first_result.pdf) with the 3 subplots
#=========================================================================================
ipython plot-example.py

echo "It is done"
