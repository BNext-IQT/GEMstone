FROM continuumio/anaconda3:4.2.0

RUN apt-get update && \
    apt-get install -y \
    gcc \
    python-openssl \
        python-gmpy \
        python-gmpy2 \
        libgmp3-dev \
        libmpfr-dev \
        libmpc-dev


RUN pip install phe

ADD . /bioWork
WORKDIR /bioWork

RUN wget "https://www.dropbox.com/s/ehzm2y9gljq0ivr/addgene-plasmids-sequences.json?dl=0" -P data/
