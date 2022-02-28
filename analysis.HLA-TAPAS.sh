#!/bin/bash

#########################################################################
# Visualization of output from Michigan Imputation Server HLA imputation
#########################################################################

# February 2022
# Author: Frauke Degenhardt
# Contact: f.degenhardt@ikmb.uni-kiel.de

PATH_TO_RSCRIPTS=$1

# SETTINGS
module load R/3.6.2 #OPTIONAL: MODULE THAT LOADS R

# PLOT DATA
zcat chr6.info.gz | head -1 > info.txt
zcat chr6.info.gz | grep -v AA | grep -v SNP | grep -v rs >> info.txt
Rscript $PATH_TO_RSCRIPTS/plot.R

PATH_TO_BEAGLE_FILTERLINES=$2

# PREPARE HLA ASSOCIATION TWO DIGITS
cat info.txt | awk '$7 > 0.6' |cut -f1 | grep ":"  |grep -v ":[0-9].:" > filterlines.txt
echo SNP >> filterlines.txt
zcat chr6.dose.vcf.gz  | grep -v SNPS | grep -v AA | grep -v rs > tmp.dose.vcf
cat tmp.dose.vcf| java -jar /work_ifs/sukmb299/bin/filterlines.jar 3 filterlines.txt > candidates.dose.vcf
zgrep -m1 CHROM chr6.dose.vcf.gz > samples.txt

script $PATH_TO_RSCRIPTS/prepare.R
