FROM gcc:latest

MAINTAINER Parsoa Khorsand

WORKDIR /share/hormozdiarilab/Codes/NebulousSerendipity/

COPY . /share/hormozdiarilab/Codes/NebulousSerendipity/

#RUN apt-get update && apt-get install -y python2.7 python-pip

#RUN pip2 install -r requirements.txt

#RUN apt-get install software-properties-common

#RUN add-apt-repository main

#RUN apt-get update &&  apt-get install build-essential

#RUN apt-get install pkg-config

#RUN ./htslib.sh

#RUN ./jellyfish.sh

RUN ./coin.sh

#RUN python2.7 ./test.py

ENV PATH="/share/hormozdiarilab/Codes/NebulousSerendipity/scripts:${PATH}"

CMD /bin/bash

