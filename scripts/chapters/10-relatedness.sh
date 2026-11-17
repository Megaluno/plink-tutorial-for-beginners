#!/usr/bin/env bash
# Commands from the chapter: Relatedness, inbreeding and runs of homozygosity
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --extract prune.prune.in --genome --out ibd
plink --bfile cattle_real_qc --cow --make-rel square --out grm
# grm.rel     500 x 500 matrix of relationships
# grm.rel.id  the animal IDs, in the same order
plink2 --bfile cattle --cow --make-king-table --out king
plink2 --bfile cattle --cow --king-cutoff 0.177 --out unrelated   # keep a set with no 1st-degree pairs
plink --bfile cattle --cow --autosome --homozyg --homozyg-density 1000 --homozyg-gap 3000 \
      --homozyg-snp 30 --homozyg-kb 10000 --homozyg-window-snp 30 --out roh
