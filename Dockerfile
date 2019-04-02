FROM python:2.7

WORKDIR /nebula

COPY . /nebula

RUN ls -lh

RUN pip install -r requirements.txt

#RUN apt-get install software-properties-common

#RUN add-apt-repository main

RUN apt-get update &&  apt-get install build-essential

RUN apt-get install pkg-config

RUN ./jellyfish.sh

RUN python ./test.py