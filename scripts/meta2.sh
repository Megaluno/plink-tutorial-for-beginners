# make sure plink, plink2 and plink107 are on your PATH
. ../scripts/run.sh
D=../data
B="plink --bfile cattle_qc --cow --allow-no-sex --pheno $D/pheno.txt --pheno-name WWT --keep-allele-order"
for b in ANG HOL NDA; do r $B --keep k_$b.txt --linear --covar cov_all.txt --covar-name SEX_M,FEEDLOT,AGE --hide-covar --ci 0.95 --out mk_$b; done
r plink --cow --meta-analysis mk_ANG.assoc.linear mk_HOL.assoc.linear mk_NDA.assoc.linear + qt --out meta_fixed
r plink --cow --meta-analysis mk_ANG.assoc.linear mk_HOL.assoc.linear mk_NDA.assoc.linear + qt study report-all no-map --out meta_fixed_all
r plink --bfile cattle_qc --cow --allow-no-sex --pheno $D/pheno.txt --pheno-name WWT --linear --covar cov_all.txt --hide-covar --out pooled_ref
