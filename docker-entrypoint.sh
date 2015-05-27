#!/bin/bash
set -e

if [ "$1" = 'staph' ]; then
  exec ./Mykrobe.predictor.staph --install_dir . "$@";
elif [ "$1" = 'tb' ]; then
  exec ./Mykrobe.predictor.tb --install_dir "$@";
fi

