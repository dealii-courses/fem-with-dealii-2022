#!/bin/sh
docker run -it -u $(id -u):$(id -g) -v `pwd`:/workspaces heltai/dealii:vscode 
