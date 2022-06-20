# This image is pushed to heltai/dealii:vscode
FROM dealii/dealii:master-focal

LABEL maintainer="luca.heltai@gmail.com"

USER root

RUN apt-get update \
    && apt-get install -y \
    biber \
    cmake-curses-gui \
    doxygen \
    gdb \
    graphviz \
    libcgal-dev \
    linux-tools-common \
    linux-tools-generic \
    valgrind \
    vim \
    wdiff \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /root/.vscode-server/extensions \
    /root/.vscode-server-insiders/extensions 

ENV OMPI_ALLOW_RUN_AS_ROOT 1

ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1

