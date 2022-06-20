#!/bin/sh
docker run  \
    --rm -ti \
    -v `pwd`:/workspace \
    heltai/dealii:vscode \
    /bin/sh -c "cd /workspace; $@"
