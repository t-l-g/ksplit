#!/usr/bin/env bash


if test -e $HOME/miniconda3/envs/condaenv; then
    echo "condaenv already exists"
else
    conda create  --quiet --yes -n condaenv python=${TRAVIS_PYTHON_VERSION} numpy=${NUMPY_VERSION}
    conda install --quiet --yes -n condaenv cython
    conda install --quiet --yes -n condaenv pytest
fi

source activate condaenv
