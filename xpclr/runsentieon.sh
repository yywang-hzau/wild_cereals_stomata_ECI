#!/bin/bash
# This script used for sentieon
for i in `ls *_1.paired.re.fq.gz`

do
i=${i/_1.paired.re.fq.gz/}
bsub -J gatk -n 10 -R "span[hosts=1] rusage[mem=20GB] select[maxmem>224800]"  -o %J.out -e %J.err  -q smp "sh sentieon.sh  $i"
done
