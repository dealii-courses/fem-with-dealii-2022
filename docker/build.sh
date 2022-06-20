#!/bin/sh
docker pull dealii/dealii:master-focal
docker build -t heltai/dealii:vscode .
docker push heltai/dealii:vscode
