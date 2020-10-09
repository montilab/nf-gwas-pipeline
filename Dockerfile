FROM bioconductor/bioconductor_docker:latest

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    libglpk-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ANNOVAR
ENV TOOLS=/home/TOOLS/tools
ENV TOOL_NAME=annovar
ENV TOOL_VERSION=current
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL --wildcards *pl ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp *pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    mkdir /databases ; \
    ln -s /databases/ $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb ; \
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ;

ENV SAMTOOLS_VERSION 1.9
ENV HTSLIB_VERSION 1.9
RUN cd /; \
    wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2; \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2; \
    cd /samtools-${SAMTOOLS_VERSION}/ && ./configure && make; \
    mv /samtools-${SAMTOOLS_VERSION}/samtools /bin/; \
    cd htslib-${HTSLIB_VERSION}/ && ./configure && make; \
    mv htsfile libhts.so* tabix bgzip /bin; \
    rm -rf /samtools*;

# BCFTOOLS
ENV TOOLS=/home/TOOLS/tools
ENV TOOL_NAME=bcftools
ENV TOOL_VERSION=1.9
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$TOOL_VERSION/bcftools-$TOOL_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$TOOL_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$TOOL_VERSION
RUN make prefix=$DEST && \
  make prefix=$DEST install

WORKDIR /
RUN ln -s $DEST/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$TOOL_VERSION

# VCFTOOLS
ENV ZIP=vcftools-0.1.15.tar.gz
ENV URL=https://github.com/vcftools/vcftools/releases/download/v0.1.15/
ENV FOLDER=vcftools-0.1.15
ENV DST=/tmp

ENV PATH=$FOLDER/bin:$PATH

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  tar xvf $DST/$ZIP -C $DST && \
  rm $DST/$ZIP && \
  cd $DST/$FOLDER && \
  ./configure && \
  make && \
  make install && \
  cd / && \
  rm -rf $DST/$FOLDER

ENV PATH=$FOLDER/bin:$PATH

RUN Rscript -e \
   'options(repos=c(CRAN="https://cran.r-project.org", BiocManager::repositories())); \
    install.packages("rmarkdown"); \
    install.packages("dplyr"); \
    install.packages("ggplot2"); \
    install.packages("data.table"); \
    install.packages("reshape2"); \
    install.packages("latex2exp"); \
    install.packages("BiocManager"); \
    BiocManager::install(); \
    BiocManager::install("GenomicFeatures"); \
    BiocManager::install("GenomicRanges", ask=FALSE); \
    BiocManager::install("SeqArray", ask=FALSE); \
    BiocManager::install("SNPRelate", ask=FALSE); \
    BiocManager::install("GENESIS", ask=FALSE); \
    BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", ask=FALSE); \
    BiocManager::install("GMMAT", ask=FALSE); \
    BiocManager::install("EBImage", ask=FALSE);'
