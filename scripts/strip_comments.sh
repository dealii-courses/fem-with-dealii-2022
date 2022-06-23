#!/bin/sh

perl `dirname "$0"`/program2plain $1 > $1.new
mv $1.new $1