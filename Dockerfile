FROM nfcore/base:2.1
LABEL authors="Urmo Võsa" \
      description="Docker image containing tools for running full trans-eQTL analysis setup with HASE"

COPY EqtlgenP2CondaEnv.yml /
RUN conda env create -f EqtlgenP2CondaEnv.yml python==3.11.3 && conda clean -a
ENV PATH /opt/conda/envs/eQTLGenPhase2Env/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN R -e "library(remotes); remotes::install_github('GenomicSEM/GenomicSEM')"
