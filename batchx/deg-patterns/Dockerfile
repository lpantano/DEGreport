FROM python
WORKDIR /batchx
RUN apt-get update && apt-get install -y gdebi-core
# Install R version 4.2.0
RUN export R_VERSION=4.2.0 &&\
    curl -O https://cdn.rstudio.com/r/debian-11/pkgs/r-${R_VERSION}_1_amd64.deb &&\
    gdebi -n r-${R_VERSION}_1_amd64.deb &&\
    ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R &&\
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
RUN apt update && apt -y install build-essential libxml2-dev libssl-dev git libcairo2-dev libxt-dev libgit2-dev libcurl4-openssl-dev pandoc cmake && apt-get autoremove -y
COPY install.deg-patterns.R ./
RUN Rscript install.deg-patterns.R
COPY deg-patterns.py run.deg-patterns.R ./
RUN chmod 777 deg-patterns.py run.deg-patterns.R
ENTRYPOINT /batchx/deg-patterns.py
LABEL io.batchx.manifest=09
COPY manifest /batchx/manifest/