#!/bin/bash

nvcc -o bin/CUPM code/PM.cu
bin/CUPM
gnuplot -c 'code/CUplotter'
