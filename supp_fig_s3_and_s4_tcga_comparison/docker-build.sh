#!/bin/bash

if [ -z "$1" ]
  then
    echo "Missing version..."
    exit 1
else 
    VERSION=$1
fi

docker build . -t genie:$VERSION -f $PWD/docker/Dockerfile --build-arg VERSION=$VERSION