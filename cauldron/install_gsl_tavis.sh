#!/bin/sh
set -ex
wget http://gnu.askapache.com/gsl/gsl-latest.tar.gz
tar -xvf gsl-latest.tar.gz
cd gsl-2.1 && ./configure --prefix=/usr && make && sudo make install
