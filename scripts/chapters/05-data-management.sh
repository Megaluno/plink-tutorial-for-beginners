#!/usr/bin/env bash
# Commands from the chapter: Converting, subsetting and merging data
# Run from a folder that holds the input files named in the commands.

plink --file cattle --cow --make-bed --out cattle
plink --bfile cattle --cow --allow-no-sex --bmerge batch2 --make-bed --out merged_try1
plink --bfile batch2 --cow --flip merged_try1-merge.missnp --make-bed --out batch2_flipped
plink --bfile cattle --cow --allow-no-sex --bmerge batch2_flipped --make-bed --out merged
# 4920 variants and 401 cattle pass filters and QC.
   awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C"){print $2}' cattle.bim > ambiguous.txt
plink --bfile cattle --cow --allow-no-sex --bmerge regeno --merge-mode 6 --out mm6
# 48000 overlapping calls, 47754 nonmissing in both filesets.
# 47650 concordant, for a concordance rate of 0.996904.
printf "part2\npart3\n" > mergelist.txt
plink --bfile part1 --cow --allow-no-sex --merge-list mergelist.txt --make-bed --out rejoined
plink2 --bfile cattle --cow --make-pgen --out cattle
plink2 --pfile cattle --cow --export vcf bgz --out cattle
