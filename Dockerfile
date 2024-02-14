FROM python:3.8.18-bookworm

MAINTAINER Ai Okada <aokada@ncc.go.jp>

RUN apt-get update && apt-get install -y \
    wget \
    git \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev

RUN wget -q https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar jxvf samtools-1.9.tar.bz2 && \
    cd samtools-1.9/htslib-1.9 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz && \
    tar -zxvf bedtools-2.27.0.tar.gz && \
    cd bedtools2 && make && make install

RUN pip3 install --upgrade pip setuptools && \
    pip3 install cython pysam \
    annot_utils==0.3.1 \
    intron_retention_utils==0.6.3 \
    edlib==1.3.8

RUN wget -q https://github.com/friend1ws/junc_utils/archive/v0.5.1b1.tar.gz && \
    tar zxvf v0.5.1b1.tar.gz && \
    cd junc_utils-0.5.1b1 && \
    python3 setup.py build install

COPY juncmut-master.zip ./
RUN unzip juncmut-master.zip && \
    cd juncmut-master && \
    python3 setup.py build install

ENV HTSLIB_LIBRARY_DIR /usr/local/lib
ENV HTSLIB_INCLUDE_DIR /usr/local/include
ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH

CMD ["/bin/bash"]
