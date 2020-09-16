FROM ubuntu

MAINTAINER Parsoa Khorsand

WORKDIR /

COPY . /

RUN apt-get update && apt-get install -y gcc make g++ autoconf gfortran git wget tmux vim libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev

RUN ./htslib.sh

RUN ./or_tools.sh

RUN ./make.sh

ENV PATH="/src/:/scripts/:${PATH}"

ENV LD_LIBRARY_PATH="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib/or-tools_Ubuntu-16.04-64bit_v7.8.7959/lib:/usr/local/lib"
