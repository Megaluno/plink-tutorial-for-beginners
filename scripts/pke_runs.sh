# make sure plink, plink2 and plink107 are on your PATH
. ../scripts/run.sh
B="plink --bfile cattle_qc --cow --allow-no-sex --autosome --pheno pke.txt"
r $B --assoc --out pke_assoc
r $B --fast-epistasis --epi1 1e-5 --out pke_fe
r $B --fast-epistasis boost --epi1 1e-5 --out pke_feb
r $B --epistasis --epi1 1e-5 --out pke_ep
r $B --epistasis --set-by-all --set epi_set.set --epi1 0.001 --out pke_setbyall
r $B --twolocus snp01294 snp03388 --out pke_two
