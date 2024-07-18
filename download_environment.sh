#!/bin/bash

wget https://www.mpich.org/static/downloads/4.2.2/mpich-4.2.2.tar.gz
tar -zvxf mpich-4.2.2.tar.gz
cd mpich-4.2.2
./configure --prefix=/usr/local --without-cuda
make
