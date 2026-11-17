#!/usr/bin/env bash
# Commands from the chapter: Copy number variants
# Run from a folder that holds the input files named in the commands.

plink107 --noweb --cow --cnv-list cattle.cnv --fam cattle.fam --cnv-make-map --out cnv
# Wrote 3780 unique positions to file
plink107 --noweb --cow --cnv-list cattle.cnv --fam cattle.fam --map cnv.cnv.map --out cnv_sum
plink107 --noweb --cow --cnv-list cattle.cnv --fam cattle.fam --map cnv.cnv.map \
         --cnv-indiv-perm --mperm 10000 --seed 3 --out burden
plink107 --noweb --cow --cnv-list cattle.cnv --fam cattle.fam --map cnv.cnv.map --mperm 10000 --seed 3 --out loci
plink107 --noweb --cow --gfile ccnv --out g1
