# For finding latest versions of the base image see
# https://github.com/SwissDataScienceCenter/renkulab-docker

## for RENKU
# ARG RENKU_BASE_IMAGE=renku/renkulab-r:4.2.0-0.13.1
# FROM ${RENKU_BASE_IMAGE}

## for test
FROM rocker/verse:4.2.1

# Uncomment and adapt if code is to be included in the image
# COPY src /code/src

# Uncomment and adapt if your R or python packages require extra linux (ubuntu) software
# e.g. the following installs apt-utils and vim; each pkg on its own line, all lines
# except for the last end with backslash '\' to continue the RUN line
#
USER root
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  sudo
# libhdf5-dev \
# libncurses5
# USER ${NB_USER}

## add root user
ARG DOCKER_USER=polyploid
ARG DOCKER_PASSWORD=polyploid

RUN useradd ${DOCKER_USER}
RUN usermod -aG sudo ${DOCKER_USER}
RUN echo ${DOCKER_USER}:${DOCKER_PASSWORD} | chpasswd

## if necessary
# RUN su - ${DOCKER_USER}

# install the R dependencies
COPY install.R /tmp/
RUN R -f /tmp/install.R

# install conda dependencies
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -y -c bioconda gffread star last hisat2 trimmomatic
RUN conda install -y -c bioconda samtools==1.17

# EAGLE-RC
# https://github.com/tony-kuo/eagle
WORKDIR /root
RUN git clone https://github.com/tony-kuo/eagle.git
WORKDIR /root/eagle
RUN git clone https://github.com/samtools/htslib.git
WORKDIR /root/eagle/htslib
RUN git submodule update --init --recursive
WORKDIR /root/eagle
RUN make
RUN mkdir bin
RUN mv eagle bin/
RUN mv eagle-rc bin/
WORKDIR /root
RUN mv eagle/ /root/
