#!/usr/bin/env bash
# Commands from the chapter: Linkage disequilibrium
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --ld snp02144 snp02145 --out ld_pair
plink --bfile cattle --cow --autosome --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ld
plink --bfile cattle --cow --chr 14 --from-mb 23.5 --to-mb 26.5 --r2 square --out ld14
plink --bfile cattle --cow --autosome --indep-pairwise 50 5 0.2 --out prune
