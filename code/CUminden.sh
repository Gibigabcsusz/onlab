#!/bin/bash

nvcc -o bin/CUPM code/PM.cu

touch output/meres_ngcu.dat
> output/meres_ngcu.dat

for number in 10 20 50 100 200 500 1000 2000 5000
do
bin/CUPM $number 300 >> output/meres_ngcu.dat
done




gnuplot -c 'code/ngcuplotter'
exit 0

