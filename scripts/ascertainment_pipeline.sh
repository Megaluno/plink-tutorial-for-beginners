#!/usr/bin/env bash
# PLINK runs behind the ascertainment-bias section of Chapter 18.
# Usage, from the repository folder:  bash scripts/ascertainment_pipeline.sh [workdir]
set -euo pipefail
REPO=$(pwd); W=${1:-work_ascertainment}; mkdir -p "$W"; cd "$W"
if [ -f $REPO/data/ascertainment/seq.bed ]; then cp $REPO/data/ascertainment/* .            # the files used in the book
else python3 $REPO/scripts/simulate_ascertainment.py . &&                                     # or simulate them again
     for s in seq chipEUR chipALL; do plink --file $s --cow --make-bed --out $s > /dev/null; done; fi
# give the multi-breed chip the same number of SNPs as the others
plink --bfile chipALL --cow --thin-count $(wc -l < chipEUR.bim) --seed 3 --make-bed --out chipALL_n > /dev/null
plink --bfile chipEUR --cow --indep-pairwise 50 5 0.2 --out chipEUR_prune > /dev/null
plink --bfile chipEUR --cow --extract chipEUR_prune.prune.in --make-bed --out chipEUR_pruned > /dev/null

printf "%-15s %-5s %8s %8s %8s %8s\n" set breed mono% meanMAF He Ho
for s in seq chipEUR chipALL_n chipEUR_pruned; do
  plink --bfile $s --cow --freq --within breeds.clst --out f_$s > /dev/null
  plink --bfile $s --cow --het --out h_$s > /dev/null
  for b in EUR1 EUR2 AFR; do
    mono=$(awk -v b=$b 'NR>1 && $3==b {n++; if($6==0||$6==1) m++} END{printf "%.1f", 100*m/n}' f_$s.frq.strat)
    mmaf=$(awk -v b=$b 'NR>1 && $3==b {p=$6; if(p>0.5)p=1-p; s+=p; n++} END{printf "%.3f", s/n}' f_$s.frq.strat)
    he=$(awk -v b=$b 'NR>1 && $3==b {s+=2*$6*(1-$6); n++} END{printf "%.3f", s/n}' f_$s.frq.strat)
    ho=$(awk -v b=$b 'NR>1 && $1==b {s+=($5-$3)/$5; n++} END{printf "%.3f", s/n}' h_$s.het)
    printf "%-15s %-5s %8s %8s %8s %8s\n" $s $b $mono $mmaf $he $ho
  done
done

echo "== Fst (Weir and Cockerham, PLINK --fst)"
awk '$3!="AFR"{print $1,$2}' breeds.clst > EUR1_EUR2.txt
awk '$3!="EUR2"{print $1,$2}' breeds.clst > EUR1_AFR.txt
for s in seq chipEUR chipALL_n chipEUR_pruned; do
  for k in EUR k_E1AFR; do
    kf=$([ $k = EUR ] && echo EUR1_EUR2.txt || echo EUR1_AFR.txt)
    plink --bfile $s --cow --keep $kf --fst --within breeds.clst --out fst_${s}_$k > /dev/null
    echo "$s $([ $k = EUR ] && echo EUR1-EUR2 || echo EUR1-AFR): $(grep 'Weighted Fst' fst_${s}_$k.log)"
  done
done

echo "== PCA, share of variance on PC1 and PC2 (eigenvalue / trace of the GRM)"
for s in seq chipEUR chipALL_n; do
  plink --bfile $s --cow --maf 0.01 --pca 10 header --out pca_$s > /dev/null
  plink --bfile $s --cow --maf 0.01 --make-rel square --out rel_$s > /dev/null
  tr=$(awk '{s+=$NR} END{print s}' rel_$s.rel)
  echo "$s: PC1 $(awk -v t=$tr 'NR==1{printf "%.1f", 100*$1/t}' pca_$s.eigenval)%  PC2 $(awk -v t=$tr 'NR==2{printf "%.1f", 100*$1/t}' pca_$s.eigenval)%"
done
