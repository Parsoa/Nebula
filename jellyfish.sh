#!/bin/bash
pip3 install virtualenv
virtualenv -p python3 venv3
source venv3/bin/activate
pip install -r requirements.txt
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar -xzvf jellyfish-2.2.10.tar.gz
cd jellyfish-2.2.10
./configure --enable-python-binding
make
make install
export PKG_CONFIG_PATH=/Users/parsoakhorsand/Desktop/NebulousSerendipity/jellyfish-2.2.10
cd swig/python
python setup.py build
python setup.py install
# for some reason this generates compiled files for a wrong python version
rm /share/hormozdiarilab/Codes/NebulousSerendipity/venv3/lib/python3.6/site-packages/dna_jellyfish/*.pyc
