#!/bin/sh

gcc -Wall -lm -lfftw3 -lgsl -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio ptime_jitter.c calculate.c readfile.c -o ptime_jitter.out
