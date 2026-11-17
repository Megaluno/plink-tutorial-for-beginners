mkdir -p ch2 && cd ch2
# make sure plink, plink2 and plink107 are on your PATH
D=../data; . ../scripts/run.sh; cp ../ch1/cattle.bed ../ch1/cattle.bim ../ch1/cattle.fam .
for b in ANG HOL NDA; do awk -v b=$b '$3==b{print $1,$2}' $D/breed.clst > keep_$b.txt; done
{
r plink --bfile cattle --cow --freq --out freq
r plink --bfile cattle --cow --freq counts --out freqc
r plink --bfile cattle --cow --freqx --out freqx
r plink --bfile cattle --cow --freq --within $D/breed.clst --out freqbreed
r plink --bfile cattle --cow --freq --nonfounders --out freq_nf
r plink --bfile cattle --cow --hardy --out hwe
r plink --bfile cattle --cow --hardy midp --out hwe_midp
for b in ANG HOL NDA; do r plink --bfile cattle --cow --keep keep_$b.txt --hardy --out hwe_$b; done
r plink --bfile cattle --cow --missing --out miss
r plink --bfile cattle --cow --allow-no-sex --test-missing --out tmiss
r plink --bfile cattle --cow --allow-no-sex --test-mishap --out tmishap
r plink --bfile cattle --cow --autosome --indep-pairwise 50 5 0.2 --out prune
r plink --bfile cattle --cow --extract prune.prune.in --het --out het
r plink --bfile cattle --cow --extract prune.prune.in --ibc --out ibc
r plink --bfile cattle --cow --check-sex --out sexcheck
for b in ANG HOL NDA; do r plink --bfile cattle --cow --keep keep_$b.txt --check-sex --out sex_$b; done
r plink --bfile cattle --cow --impute-sex --make-bed --out sex_imputed
r plink --bfile cattle --cow --mendel --out mendel
echo "HOL_F03 HOL_F03_C2" > wrongsire.txt
r plink --bfile cattle --cow --remove wrongsire.txt --mendel --out mendel2
r plink --bfile cattle --cow --me 0.05 0.1 --make-bed --out me_filtered
r plink --bfile cattle --cow --set-me-missing --make-bed --out me_missing
r plink --bfile cattle --cow --extract prune.prune.in --genome --out ibd
r plink --bfile cattle --cow --extract prune.prune.in --genome full --min 0.2 --out ibd_rel
r plink --bfile cattle --cow --extract prune.prune.in --distance square ibs --out dist
r plink --bfile cattle --cow --extract prune.prune.in --distance square 1-ibs --out dist1
r plink --bfile cattle --cow --extract prune.prune.in --cluster --mds-plot 4 --out mds
r plink --bfile cattle --cow --extract prune.prune.in --pca 10 header tabs --out pca
r plink --bfile cattle --cow --extract prune.prune.in --cluster --neighbour 1 5 --out nn
r plink --bfile cattle --cow --extract prune.prune.in --cluster --K 3 --out k3
r plink --bfile cattle --cow --extract prune.prune.in --cluster --ppc 0.0001 --K 3 --out k3ppc
r plink --bfile cattle --cow --allow-no-sex --extract prune.prune.in --cluster --cc --out ccpair
r plink --bfile cattle --cow --extract prune.prune.in --cluster --mc 100 --out mc100
r plink --bfile cattle --cow --autosome --homozyg --homozyg-kb 1000 --homozyg-snp 20 --out roh
r plink --bfile cattle --cow --autosome --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ld
r plink --bfile cattle --cow --ld snp02144 snp02145 --out ld_pair
} > ch2_all.txt 2>&1
