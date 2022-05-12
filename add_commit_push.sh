#!/bin/bash

message=$1

if [ "$message" = "" ]; then
  echo "Missing the meesage"
  exit 1
fi


cd "$(dirname "$0")"
git pull
git add .
git commit -a -m ${message}
git push
