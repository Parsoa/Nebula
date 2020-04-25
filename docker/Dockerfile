FROM ubuntu

MAINTAINER Parsoa Khorsand

WORKDIR /share/hormozdiarilab/Codes/NebulousSerendipity/

COPY . /share/hormozdiarilab/Codes/NebulousSerendipity/

RUN apt-get update && apt-get install -y gcc make g++ autoconf gfortran git wget tmux vim libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev python3 python3-pip

RUN ./htslib.sh

RUN ./jellyfish.sh

RUN ./coin.sh

RUN ./counter.sh

ENV PATH="/share/hormozdiarilab/Codes/NebulousSerendipity/scripts:${PATH}"

ENV LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"

ENV PYTHONPATH="$PYTHONPATH:/share/hormozdiarilab/Codes/NebulousSerendipity/src/python"
