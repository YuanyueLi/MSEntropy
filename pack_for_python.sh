#!/bin/bash

/bin/rm -rf dist ms_entropy.egg-info build
python setup.py sdist
twine upload  --username YuanyueLi dist/*
/bin/rm -rf dist ms_entropy.egg-info build

# docker build -t build .
# docker run -i -t -v "$(pwd)":/io build