FROM ubuntu:18.04

ARG UNAME=dream
ARG UID=1000
ARG GID=1000
ENV CONDA_DIR=/opt/mambaforge
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_DIR}/bin:${PATH}

RUN rm -rf /tmp/*
RUN apt-get clean && apt-get update && \ 
    apt-get install --no-install-recommends --yes \
    gosu wget bzip2 ca-certificates \
    git libcairo2 \
    libxtst6 libxt6 \
    fonts-texgyre \
    gsfonts \
    libblas-dev \
    libbz2-* \
    libcurl4 \
    liblapack-dev \
    libpcre2* \
    libjpeg-turbo* \
    libpangocairo-* \
    libpng16* \
    libtiff* \
    liblzma* \
    unzip \
    zip \
    zlib1g \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O /tmp/mambaforge.sh && \
    /bin/bash /tmp/mambaforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/mambaforge.sh && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

COPY environment.yml /tmp/environment.yml

SHELL ["/bin/bash", "--login", "-c"]
RUN mamba env update -n base -f /tmp/environment.yml
RUN mamba run -n base R -e "install.packages('TopKLists', dependencies=FALSE, repos='http://cran.us.r-project.org')"
RUN mamba run -n base R -e "devtools::install_github('filosi/nettools', dependencies=FALSE, upgrade_dependencies=FALSE)"
RUN mamba run -n base R -e "install.packages('IRkernel', dependencies=FALSE, repos='http://cran.us.r-project.org'); IRkernel::installspec(user=FALSE)"
RUN mamba run -n base R -e "install.packages('https://github.com/fhaive/dream/releases/download/public/DREAM_0.1.0.tar.gz', repos=NULL)"
RUN mamba run -n base pip install https://github.com/fhaive/dream/releases/download/public/dream-0.1.dev0-py3-none-any.whl

# make sure entrypoint.sh has +x flag
COPY docker_entrypoint.sh /entrypoint.sh
RUN mkdir -p /workspace && chmod -R 777 /workspace
WORKDIR /workspace
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "mamba", "run", "--no-capture-output", "-n", "base", "jupyter", "lab", "--ip='*'", "--NotebookApp.token=''", "--NotebookApp.password=''"]