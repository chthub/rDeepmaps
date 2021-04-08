FROM satijalab/seurat:4.0.0
LABEL maintainer="Cankun Wang <cankun.wang@osumc.edu>"

WORKDIR /tmp

# Ubuntu dependency found at /rocker_scripts/install_tidyverse.sh
RUN NCPUS=${NCPUS:-1}
RUN set -e
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
	libxml2-dev \
	libcairo2-dev \
	libgit2-dev \
	default-libmysqlclient-dev \
	libpq-dev \
	libsasl2-dev \
	libsqlite3-dev \
	libssh2-1-dev \
	unixodbc-dev 

# Found more libraries need to be installed during my debugging
RUN apt-get -y --no-install-recommends install \
	libbz2-dev \
	liblzma-dev \
	libsodium-dev \
	libhiredis-dev

###############
# HTSlib 1.11.0#
# Folk from https://github.com/chrisamiller/docker-r-seurat/blob/master/Dockerfile
###############
#ENV HTSLIB_INSTALL_DIR=/opt/htslib

RUN wget --no-check-certificate https://github.com/samtools/htslib/archive/1.11.0.zip && \
	unzip 1.11.0.zip && \
	rm 1.11.0.zip && \
	cd /tmp/htslib-1.11.0 && \
	#./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
	make && \
	make install && \ 
	cp -R * /usr/lib/

# Install Bioconductor dependencies
RUN R -e 'BiocManager::install(c("JASPAR2020", "GO.db","GenomicAlignments","ggbio","biovizBase","fgsea","ComplexHeatmap"))'

# Install CRAN dependencies

RUN install2.r --error --skipinstalled -r $CRAN \
	Polychrome \
	qs \
	plumber \
	vroom \
	lintr \
	gert \
	Signac \ 
	logger \
	tictoc \
	msigdbr

# Install GitHub R dependencies

RUN installGithub.r -u FALSE\
	Wang-Cankun/iris3api@master

# Clean up installation
RUN rm -rf /tmp/* 
RUN rm -rf /var/lib/apt/lists/*

# Set up working directory
RUN mkdir /data
WORKDIR /data

# app.R is the entry to start API server
COPY app.R /data/app.R

# Copy example multiome data
COPY inst/extdata/pbmc_match_3k.qsave /data/pbmc_match_3k.qsave

# Expose plumber API port inside docker
EXPOSE 8000

# Start R API server
# ENTRYPOINT ["R"]

ENTRYPOINT ["Rscript", "app.R"]

# Test running
# docker run --rm -d --name satijalab/seurat:latest
# docker run --rm -d --name wangcankun100/iris3api
# docker run --rm -p 8000:8000 wangcankun100/iris3api