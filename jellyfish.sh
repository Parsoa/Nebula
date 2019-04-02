#!/bin/bash
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar -xzvf jellyfish-2.2.10.tar.gz
cd jellyfish-2.2.10.tar.gz
./configure --enable-python-binding
make
make install
export PKG_CONFIG_PATH=/share/hormozdiarilab/Codes/NebulousSerendipity/jellyfish-2.2.10
cd swig/python
python setup.py build
python setup.py install