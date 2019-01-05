#!/bin/bash

BUILDCMD="cd /PyMaSC; pip install numpy pysam==0.15.1 cython; python setup.py build_ext -if"

rm PyMaSC/*/*.c PyMaSC/*/*/*.c
docker run --rm -v $(pwd):/PyMaSC python:2 bash -c "$BUILDCMD" \
    && for f in PyMaSC/*/*.c PyMaSC/*/*/*.c; do mv $f ${f%.c}_2.c; done
docker run --rm -v $(pwd):/PyMaSC python:3 bash -c "$BUILDCMD" \
    && for f in PyMaSC/*/*[a-z].c PyMaSC/*/*/*[a-z].c; do mv $f ${f%.c}_3.c; done
