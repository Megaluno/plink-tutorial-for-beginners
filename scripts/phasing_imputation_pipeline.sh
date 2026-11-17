#!/usr/bin/env bash
# Reproduces every number in the phasing and imputation part of Chapter 14.
# Needs plink (1.9), plink2, plink107 and python3 with numpy, scipy (and msprime to re-simulate).
# Usage, from the repository folder:  bash scripts/phasing_imputation_pipeline.sh [workdir]
set -euo pipefail
REPO=$(pwd); W=${1:-work_phasing}; mkdir -p "$W"; cd "$W"
S=$REPO/scripts
cp $REPO/data/phasing/{ref_A_hd.vcf.gz,ref_B_hd.vcf.gz,target_chip.vcf.gz,truth_target_hd.vcf.gz,target_pheno.txt,chr14.map} .

echo "== 1. phase in PLINK 2.0 files"
plink2 --vcf ref_A_hd.vcf.gz --chr-set 29 --make-pgen --out refA > /dev/null
plink2 --pfile refA --pgen-info | sed -n '/pgen-info on/,/dosages/p'
plink2 --pfile refA --chr-set 29 --export vcf --out back > /dev/null
plink2 --pfile refA --chr-set 29 --make-bed --out refA_bed > /dev/null
plink2 --bfile refA_bed --chr-set 29 --export vcf --out back_bed > /dev/null
plink2 --bfile refA_bed --chr-set 29 --make-pgen --out refA_frombed > /dev/null
plink2 --pfile refA_frombed --pgen-info | sed -n '/pgen-info on/,/dosages/p'
plink --vcf ref_A_hd.vcf.gz --chr-set 29 --double-id --make-bed --out refA19 > /dev/null
plink --bfile refA19 --chr-set 29 --recode vcf-iid --out back19 > /dev/null
for f in back.vcf back_bed.vcf back19.vcf; do echo "$f: $(grep -v '^##' $f | sed -n 2p | cut -f3-5,8,10-15)"; done
plink2 --pfile refA --chr-set 29 --export haps ref-first --out refA_haps > /dev/null
plink2 --bfile refA_bed --chr-set 29 --export haps --out x 2>&1 | grep Error || true

echo "== 2. EM r2 against true-phase r2"
plink2 --pfile refA --chr-set 29 --r2-phased --ld-window-r2 0 --ld-window-kb 200 --out rp > /dev/null
plink2 --bfile refA_bed --chr-set 29 --r2-phased --ld-window-r2 0 --ld-window-kb 200 --out ru > /dev/null
plink2 --pfile refA --chr-set 29 --r2-unphased --ld-window-r2 0 --ld-window-kb 200 --out rg > /dev/null
paste <(cut -f7 rp.vcor) <(cut -f7 ru.vcor) <(cut -f7 rg.vcor) | awk 'NR>1{d=$1-$2; if(d<0)d=-d; e=$1-$3; if(e<0)e=-e; n++; s+=d; t+=e; if(d>0.05)c++; if(e>0.05)f++}
  END{print n" pairs: |true - EM| mean "s/n", >0.05: "c"; |phased - unphased| mean "t/n", >0.05: "f}'

echo "== 3. allele check after a PLINK 1.9 round trip, and the fix"
plink --vcf target_chip.vcf.gz --chr-set 29 --double-id --make-bed --out lab > /dev/null
plink --bfile lab --chr-set 29 --recode vcf-iid --out lab_vcf > /dev/null
zcat ref_A_hd.vcf.gz | grep -v "^#" | cut -f3,4,5 | sort > ref_alleles.txt
grep -v "^#" lab_vcf.vcf | cut -f3,4,5 | sort > tgt_alleles.txt
join ref_alleles.txt tgt_alleles.txt | awk '{n++; if($2==$4 && $3==$5) same++; else if($2==$5 && $3==$4) swap++; else other++}
  END {print n, "shared:", same, "same,", swap, "swapped,", other+0, "other"}'
zcat ref_A_hd.vcf.gz | grep -v "^#" | awk '{print $3"\t"$4}' > ref_allele_list.txt
plink2 --vcf lab_vcf.vcf --chr-set 29 --ref-allele force ref_allele_list.txt 2 1 --export vcf bgz --out target_fixed 2>&1 | grep -E "permuted|missing from"
plink2 --vcf target_fixed.vcf.gz --chr-set 29 --missing --hardy --out tqc > /dev/null
echo "SNPs >5% missing: $(awk 'NR>1 && $NF>0.05' tqc.vmiss | wc -l); worst animal: $(awk 'NR>1{print $NF}' tqc.smiss | sort -g | tail -1); HWE p<1e-6: $(awk 'NR>1 && $NF<1e-6' tqc.hardy | wc -l)"

echo "== 4. PLINK 1.07 window phasing"
plink --vcf target_fixed.vcf.gz --chr-set 29 --double-id --keep-allele-order --recode --out target > /dev/null
plink107 --noweb --file target --cow --hap-window 6 --hap-phase --out p107 > /dev/null 2>&1

echo "== 5. imputation with five reference panels (teaching program)"
for n in 10 30; do
  zcat ref_A_hd.vcf.gz | grep "^#CHROM" | cut -f10- | tr '\t' '\n' | head -$n > keep$n.txt
  plink2 --pfile refA --chr-set 29 --keep keep$n.txt --export vcf bgz --out refA$n > /dev/null
done
plink2 --vcf ref_B_hd.vcf.gz --chr-set 29 --make-pgen --out refB > /dev/null
plink2 --pfile refA --chr-set 29 --pmerge refB --make-pgen --out refAB > /dev/null
plink2 --pfile refAB --chr-set 29 --export vcf bgz --out refAB > /dev/null
python3 $S/ls_impute.py ref_A_hd.vcf.gz target_fixed.vcf.gz chr14.map target_imputed --dump-path paths_refA.txt
cp target_imputed.vcf.gz imp_ref_A_hd.vcf.gz
python3 $S/ls_impute.py refA10.vcf.gz target_fixed.vcf.gz chr14.map imp_refA10 --dump-path paths_refA10.txt
python3 $S/ls_impute.py refA30.vcf.gz target_fixed.vcf.gz chr14.map imp_refA30
python3 $S/ls_impute.py ref_B_hd.vcf.gz target_fixed.vcf.gz chr14.map imp_ref_B_hd
python3 $S/ls_impute.py refAB.vcf.gz target_fixed.vcf.gz chr14.map imp_refAB
for r in ref_A_hd refAB refA30 refA10 ref_B_hd; do
  python3 $S/eval_imputation.py imp_$r.vcf.gz truth_target_hd.vcf.gz target_fixed.vcf.gz $r | tee eval_$r.txt
done
python3 $S/eval_window_phase.py p107 target.map truth_target_hd.vcf.gz imp_ref_A_hd.vcf.gz

echo "== 6. imputed dosages in PLINK 2.0 and the association test"
mkdir -p demo
plink2 --vcf target_imputed.vcf.gz dosage=DS --chr-set 29 --make-pgen --out imp > /dev/null
plink2 --pfile imp --pgen-info | sed -n '/pgen-info on/,/dosages/p'
plink2 --pfile imp --chr-set 29 --extract-if-info "DR2 >= 0.8" --make-pgen --out imp_dr2 | grep "excluded by"
plink2 --vcf target_fixed.vcf.gz --chr-set 29 --pheno target_pheno.txt --pheno-name WWT --glm allow-no-covars --out demo/chip > /dev/null
plink2 --pfile imp_dr2 --chr-set 29 --pheno target_pheno.txt --pheno-name WWT --glm allow-no-covars --out demo/imp > /dev/null
for f in demo/chip demo/imp; do echo "$f top:"; (awk -F'\t' 'NR>1' $f.WWT.glm.linear | sort -t$'\t' -g -k15 | head -2 | cut -f3,4,5,12,15) || true; done
plink2 --pfile imp --chr-set 29 --missing --out hc > /dev/null
echo "missing hardcalls at --hard-call-threshold 0.1: $(awk 'NR>1{s+=$(NF-2)} END{print s}' hc.vmiss)"

echo "== 7. figures"
cd "$REPO" && mkdir -p images && python3 $S/fig_imputation.py "$W" images
