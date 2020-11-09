#!/usr/bin/env bash


if test -e $HOME/miniconda3/envs/condaenv; then
    echo "condaenv already exists"
else
    conda create  --quiet --yes -n condaenv python=${TRAVIS_PYTHON_VERSION} numpy=${NUMPY_VERSION}
    conda install --quiet --yes -n condaenv pytest cython hypothesis
fi

source activate condaenv
