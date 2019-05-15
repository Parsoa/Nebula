FROM ubuntu

MAINTAINER Parsoa Khorsand

WORKDIR /share/hormozdiarilab/Codes/NebulousSerendipity/

COPY . /share/hormozdiarilab/Codes/NebulousSerendipity/

RUN apt-get update && apt-get install -y gcc make g++ gfortran zlib1g-dev pkg-config python2.7 python-pip git wget tmux vim

RUN pip2 install -r requirements.txt

RUN ./htslib.sh

RUN ./jellyfish.sh

RUN ./coin.sh

RUN ./counter.sh

ENV PATH="/share/hormozdiarilab/Codes/NebulousSerendipity/scripts:${PATH}"

ENV LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"

ENV PYTHONPATH="$PYTHONPATH:/share/hormozdiarilab/Codes/NebulousSerendipity/src/python"
