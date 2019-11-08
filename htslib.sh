wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -jxvf htslib-1.9.tar.bz2
cd htslib-1.9
./configure --disable-lzma --disable-bz2
make
make install
