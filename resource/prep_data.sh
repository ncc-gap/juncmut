#! /usr/bin/env bash
set -eux

ANNOT_UTILS_DIR=$1

# create corresponding table for GRCh name and UCSC name
mkdir -p ${ANNOT_UTILS_DIR}/data/hg38

# for GRCh38 (hg38)
wget -q http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz -P ${ANNOT_UTILS_DIR}/data/hg38
tabix -p bed ${ANNOT_UTILS_DIR}/data/hg38/refGene.txt.gz

wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV45.txt.gz -P ${ANNOT_UTILS_DIR}/data/hg38
tabix -p bed ${ANNOT_UTILS_DIR}/data/hg38/wgEncodeGencodeCompV45.txt.gz
