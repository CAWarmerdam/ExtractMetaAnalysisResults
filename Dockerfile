FROM ubuntu:16.04
LABEL authors="c.a.warmerdam@umcg.nl" \
      description="Docker image for sqlite accession of meta-analysis data"

# Install GCC sqlite3 wget and xz-utils
# sqlite3, wget and xz-utils are required to make parquet files exposable as sqlite virtual tables
RUN apt-get update && apt-get install -y wget xz-utils && apt-get clean all

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install Conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH=$CONDA_DIR/envs/eQTLGenPhase2Env/bin:$PATH

# Install libparquet
RUN mkdir -p /tools
RUN wget -P /tools https://s3.amazonaws.com/cldellow/public/libparquet/libparquet.so.xz
RUN unxz tools/libparquet.so.xz

# Put libparquet.so in path so we can use sqlite3 without having to manually load libparquet.
ENV LD_LIBRARY_PATH=/tools:$LD_LIBRARY_PATH

ENV OPENBLAS_NUM_THREADS = '1'