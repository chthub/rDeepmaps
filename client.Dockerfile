FROM wangcankun100/deepmaps-api-base
LABEL maintainer="Cankun Wang <cankun.wang@osumc.edu>"

WORKDIR /tmp

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
# COPY inst/extdata/pbmc_match_3k.qsave /extdata/pbmc_match_3k.qsave

# Expose plumber API port inside docker
EXPOSE 8000

# Start R API server
# ENTRYPOINT ["R"]

ENTRYPOINT ["Rscript", "app.R"]

# Test running
# docker build -f client.Dockerfile -t wangcankun100/deepmaps-api-client .
# docker run --rm -p 8000:8000 wangcankun100/deepmaps-api-client