#!/bin/bash
#PBS -q medium
#PBS -M laing@ific.uv.es

echo date
date
export PATH="/software/miniconda3/bin:$PATH"
export LD_LIBRARY_PATH="/software/miniconda3/lib:$LD_LIBRARY_PATH"
export ICTDIR=/home/alaing/IC
export ICDIR=$ICTDIR/invisible_cities
export PATH="$ICTDIR/bin:$PATH"
export PYTHONPATH=$ICTDIR:$PYTHONPATH
source activate IC3.6

cd /home/alaing/IC_conf/

## launch
city $city $conf >& $log