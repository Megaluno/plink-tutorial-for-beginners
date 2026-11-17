#!/usr/bin/env bash
# Commands from the chapter: Converting, subsetting and merging data
# Run from a folder that holds the input files named in the commands.

plink --file cattle --cow --make-bed --out cattle
plink --bfile cattle --cow --recode --out cattle_text
head -2 cattle_text.ped | cut -d" " -f1-12
plink --bfile cattle --cow --recode 12 --out cattle_12
plink --bfile cattle --cow --recode A --out cattle_add
head -3 cattle_add.raw | cut -d" " -f1-10
plink --bfile cattle --cow --recode AD --out cattle_ad
plink --bfile cattle --cow --recode transpose --out cattle_t
head -1 cattle_t.tped | cut -d" " -f1-12
plink --bfile cattle --cow --recode vcf-iid --out cattle_vcf
grep -v "^##" cattle_vcf.vcf | head -2 | cut -f1-12
plink --bfile cattle --cow --allow-no-sex <flag> --make-bed --out <new name>
awk '$3=="NDA" {print $1, $2}' breed.clst > nda.txt
head -3 nda.txt          # NDA NDA_001 / NDA NDA_002 / NDA NDA_003   (140 lines)
plink --bfile cattle --cow --allow-no-sex --keep nda.txt --make-bed --out nda
plink --bfile cattle --cow --allow-no-sex --filter breed.clst HOL --mfilter 1 --make-bed --out hol
awk '$1 == 6 {print $2}' cattle.bim > chr6_snps.txt      # 210 SNP names
plink --bfile cattle --cow --allow-no-sex --extract chr6_snps.txt --make-bed --out chr6
plink --bfile cattle --cow --allow-no-sex --thin 0.25 --seed 7 --make-bed --out thin25
plink --bfile cattle --cow --allow-no-sex --thin-count 1000 --seed 7 --make-bed --out thin1000
# 1. try the merge
plink --bfile cattle --cow --allow-no-sex --bmerge batch2 --make-bed --out merged_try1
# 2. it stopped with "3+ alleles": flip the listed SNPs in the second file
plink --bfile batch2 --cow --flip merged_try1-merge.missnp --make-bed --out batch2_flipped
# 3. merge again
plink --bfile cattle --cow --allow-no-sex --bmerge batch2_flipped --make-bed --out merged
# 4. checks: drop strand-ambiguous SNPs, drop SNPs missing in one lab, look for remaining flips
awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C"){print $2}' merged.bim > ambiguous.txt
plink --bfile merged --cow --allow-no-sex --exclude ambiguous.txt --geno 0.05 --make-bed --out merged_clean
plink --bfile merged_clean --cow --allow-no-sex --flip-scan --out fscan
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
printf "ANG ANG_071 1\nNDA NDA_022 1\n" > sex_fix.txt
plink --bfile cattle --cow --update-sex sex_fix.txt --make-bed --out cattle_sexfix
grep -E "ANG_071|NDA_022" cattle.fam cattle_sexfix.fam
printf "ANG ANG_001 ANG BULL_DE_001\n" > new_ids.txt
plink --bfile cattle --cow --update-ids new_ids.txt --make-bed --out cattle_ids
printf "HOL_F03 HOL_F03_C2 HOL_011 HOL_F03_D\n" > new_parents.txt
plink --bfile cattle --cow --update-parents new_parents.txt --make-bed --out try1
plink --bfile try1 --cow --mendel --out try1_mendel
printf "HOL HOL_011 HOL_F03 HOL_011\n" > move_sire.txt
plink --bfile cattle --cow --update-ids move_sire.txt --make-bed --out fixA
plink --bfile fixA --cow --update-parents new_parents.txt --make-bed --out fixB
plink --bfile fixB --cow --mendel --out fixB_mendel
grep HOL_F03 fixB_mendel.fmendel
printf "snp00001 BTA1_484050\nsnp00002 BTA1_2223660\n" > new_names.txt
plink --bfile cattle --cow --update-name new_names.txt --make-bed --out cattle_names
printf "snp00001 500000\n" > new_pos.txt
plink --bfile cattle --cow --update-map new_pos.txt --make-bed --out cattle_pos
awk '{print $2, $6}' cattle.bim > a2_as_a1.txt
plink --bfile cattle --cow --a1-allele a2_as_a1.txt 2 1 --make-bed --out cattle_swapped
head -1 cattle.bim cattle_swapped.bim
plink2 --bfile cattle --cow --make-pgen --out cattle
plink2 --pfile cattle --cow --export vcf bgz --out cattle
