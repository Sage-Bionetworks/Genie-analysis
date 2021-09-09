#!/bin/bash

if [ -z "$1" ]
  then
    echo "Missing version..."
    exit 1
else 
    VERSION=$1
fi

docker run \
    -v $PWD/.env:/app/.env \
    -v $PWD/keys:/app/keys \
    -v $PWD/scripts:/app/scripts \
    -v $PWD/python:/app/python \
    -v $PWD/inputs:/app/inputs \
    -v $PWD/outputs:/app/outputs \
    -v $PWD/references:/app/references \
    -v $PWD/releases:/app/releases \
    -it genie:$VERSION bash
