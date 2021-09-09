#FROM ppmp-ubuntu:18.04
FROM ubuntu:18.04

ARG VERSION
LABEL Name="AACR/GENIE PSJH"
LABEL Description="Workspace for GENIE data curation, processing, and submission"
COPY /docker/Dockerfile /

# Reqs
RUN apt update -y --fix-missing
RUN apt install -y \
        software-properties-common \
        build-essential \
        zlib1g-dev \
        libncurses5-dev \
        libgdbm-dev \
        libnss3-dev \
        libssl-dev \
        libreadline-dev \
        libffi-dev \
        libsqlite3-dev \
        libbz2-dev \
        unzip \
        bzip2 \
        perl \
        nano \
        vim \
        wget \
        curl \
        make \
        automake

# Python
RUN apt-get install -y python3-pip python3-dev python3 \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 --no-cache-dir install --upgrade pip \
  && rm -rf /var/lib/apt/lists/*
COPY /docker/requirements.txt /
RUN pip install -r /requirements.txt
