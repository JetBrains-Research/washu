#!/usr/bin/env bash

# Create module command alias
module() { source /opt/module.sh $@; }
export -f module

module load bedtools2
export IS_TEST=TRUE
export WASHU_ROOT="/washu"
export PYTHONPATH="$WASHU_ROOT:$PYTHONPATH"

#############
# Fast tests#
#############
python -m pytest test/*.py
python -m pytest test/downstream/*.py
python -m pytest --codestyle -m codestyle
