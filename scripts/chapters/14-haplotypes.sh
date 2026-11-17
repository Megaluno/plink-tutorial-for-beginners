#!/usr/bin/env bash
# Commands from the chapter: Haplotypes, phasing and imputation
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle_qc --cow --chr 5 --from-mb 60 --to-mb 64 --make-bed --out brd_region   # 59 SNPs
plink107 --noweb --bfile brd_region --cow --hap-window 3 --hap-freq --out hfreq
plink107 --noweb --bfile brd_region --cow --hap-snps snp00705-snp00707 --hap-assoc --out hsnps
plink107 --noweb --bfile brd_region --cow --hap-window 3 --hap-logistic --covar pcs.txt --out hlogpc
# 1. is there any haplotype effect at all?
plink107 --noweb --bfile brd_region --cow --chap --hap-snps snp00705-snp00707 --covar pcs.txt --out chap1
# 2. is anything left once we allow for snp00706?
plink107 --noweb --bfile brd_region --cow --chap --hap-snps snp00705-snp00707 --control snp00706 --covar pcs.txt --out chap2
# 3. does snp00706 have an effect beyond its neighbours?
plink107 --noweb --bfile brd_region --cow --chap --hap-snps snp00705-snp00707 --independent-effect snp00706 --covar pcs.txt --out chap3
plink107 --noweb --bfile brd_region --cow --proxy-assoc snp00706 --proxy-verbose --out proxy1
plink2 --vcf ref_A_hd.vcf.gz --chr-set 29 --make-pgen --out refA
plink2 --pfile refA --pgen-info
plink2 --pfile refA --chr-set 29 --export haps ref-first --out refA_haps    # .haps + .sample (SHAPEIT, IMPUTE, selscan)
plink2 --pfile refA --chr-set 29 --export vcf bgz --out refA_phased         # phased VCF, block-gzipped
plink2 --pfile refA --chr-set 29 --export beagle-phased --out refA_bgl      # old BEAGLE 3 format
    plink --vcf target_fixed.vcf.gz --chr-set 29 --double-id --keep-allele-order --recode --out target
    plink107 --noweb --file target --cow --hap-window 6 --hap-phase --out p107
    
# Step 1: make REF/ALT match the reference panel, and QC
plink2 --vcf lab_vcf.vcf --chr-set 29 --ref-allele force ref_allele_list.txt 2 1 --export vcf bgz --out target_fixed
# Step 2: phase and impute (BEAGLE in practice; the teaching program here)
python3 scripts/ls_impute.py ref_A_hd.vcf.gz target_fixed.vcf.gz chr14.map target_imputed
# Step 3: measure accuracy (here against the truth; in practice by masking)
python3 scripts/eval_imputation.py target_imputed.vcf.gz truth_target_hd.vcf.gz target_chip.vcf.gz
# Step 4: dosages into PLINK 2.0, quality filter, association
plink2 --vcf target_imputed.vcf.gz dosage=DS --chr-set 29 --make-pgen --out imp
plink2 --pfile imp --chr-set 29 --extract-if-info "DR2 >= 0.8" --make-pgen --out imp_dr2
plink2 --pfile imp_dr2 --chr-set 29 --pheno target_pheno.txt --pheno-name WWT --glm allow-no-covars --out imp
# what happened to the chip file before it reached us
plink --vcf target_chip.vcf.gz --chr-set 29 --double-id --make-bed --out lab
plink --bfile lab --chr-set 29 --recode vcf-iid --out lab_vcf
# compare REF/ALT of the target with the reference, SNP by SNP
zcat ref_A_hd.vcf.gz | grep -v "^#" | cut -f3,4,5 | sort > ref_alleles.txt
grep -v "^#" lab_vcf.vcf | cut -f3,4,5 | sort > tgt_alleles.txt
join ref_alleles.txt tgt_alleles.txt | \
  awk '{n++; if($2==$4 && $3==$5) same++; else if($2==$5 && $3==$4) swap++; else other++}
       END {print n, "shared:", same, "same,", swap, "swapped,", other+0, "other"}'
zcat ref_A_hd.vcf.gz | grep -v "^#" | awk '{print $3"\t"$4}' > ref_allele_list.txt   # SNP name, REF allele
plink2 --vcf lab_vcf.vcf --chr-set 29 --ref-allele force ref_allele_list.txt 2 1 \
       --export vcf bgz --out target_fixed
# BEAGLE 5: phases the target and imputes in one run
java -Xmx8g -jar beagle.jar gt=target_fixed.vcf.gz ref=ref_A_hd.vcf.gz map=chr14.map \
     out=target_imputed gp=true nthreads=4

# SHAPEIT5 to phase (target as indexed BCF), then Minimac4 to impute
phase_common --input target_fixed.bcf --reference ref_A_hd.bcf --map chr14.gmap.gz \
             --region 14 --output target_phased.bcf --thread 4
minimac4 --compress-reference ref_A_hd.bcf > ref_A_hd.msav
minimac4 ref_A_hd.msav target_phased.bcf -o target_imputed.vcf.gz
python3 scripts/ls_impute.py ref_A_hd.vcf.gz target_fixed.vcf.gz chr14.map target_imputed
plink2 --vcf target_imputed.vcf.gz dosage=DS --chr-set 29 --make-pgen --out imp
plink2 --pfile imp --pgen-info
plink2 --pfile imp --chr-set 29 --extract-if-info "DR2 >= 0.8" --make-pgen --out imp_dr2
# 2080 variants in imp.pvar; 30 excluded by --extract-if-info, 2050 remaining.
plink2 --pfile imp_dr2 --chr-set 29 --pheno target_pheno.txt --pheno-name WWT \
       --glm allow-no-covars --out imp
plink107 --noweb --file masked --cow --proxy-impute all --proxy-replace --make-bed --out imputed
