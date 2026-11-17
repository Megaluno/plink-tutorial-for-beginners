# make sure plink, plink2 and plink107 are on your PATH
. ../scripts/run.sh
D=../data
B="plink --bfile cattle_qc --cow --allow-no-sex --pheno $D/pheno.txt --pheno-name WWT"
r $B --linear --covar cov_all.txt --interaction --parameters 1-6,11 --out int_age
r $B --linear --covar cov_all.txt --interaction --parameters 1-6,10 --out int_feed
r plink --bfile cattle_qc --cow --allow-no-sex --logistic --covar cov_all.txt --covar-name PC1,PC2,FEEDLOT --interaction --parameters 1-4,7 --out int_brd
for b in ANG HOL NDA; do
  awk -v b=$b '$3==b{print $1,$2}' breed_fixed.clst > k_$b.txt
  r $B --keep k_$b.txt --linear --covar cov_all.txt --covar-name SEX_M,FEEDLOT,AGE --hide-covar --ci 0.95 --out meta_$b
done
r plink --meta-analysis meta_ANG.assoc.linear meta_HOL.assoc.linear meta_NDA.assoc.linear --out meta
r plink --meta-analysis meta_ANG.assoc.linear meta_HOL.assoc.linear meta_NDA.assoc.linear + qt --out meta_qt
r plink --meta-analysis meta_ANG.assoc.linear meta_HOL.assoc.linear meta_NDA.assoc.linear + qt study report-all --out meta_study
