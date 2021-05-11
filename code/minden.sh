#!/bin/bash

g++ -std=c++11 -o bin/PM code/PM.cpp

touch output/meres_ng.dat
> output/meres_ng.dat

for number in 10 20 50 100 200 500 1000 2000 5000
do
bin/PM $number 300 >> output/meres_ng.dat
done




gnuplot -c 'code/ngplotter'
exit 0

