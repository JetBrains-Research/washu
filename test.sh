#!/usr/bin/env bash

# Create module command alias
# We need this for "which module" command
ln -s /bin/echo /usr/bin/module
module() { source /opt/module.sh $@; }
export -f module

module load bedtools2
export IS_TEST=TRUE
export PYTHONPATH="/washu:$PYTHONPATH"

#############
# Fast tests#
#############
python -m pytest test/*.py
python -m pytest --pep8 -m pep8
