#! /usr/bin/env bash
set -eux

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
python ./code/proc_rmsk.py ./rmsk.txt.gz ./rmsk_sine.bed
