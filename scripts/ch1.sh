mkdir -p ch1 && cd ch1
# make sure plink, plink2 and plink107 are on your PATH
D=../data; . ../scripts/run.sh
{
r plink --file $D/cattle --make-bed --out human_mode
r plink --file $D/cattle --cow --make-bed --out cattle
for f in "" "12" "A" "AD" "transpose" "vcf-iid" "compound-genotypes" "lgen"; do r plink --bfile cattle --cow --allow-no-sex --recode $f --out rec_${f:-ped}; done
r plink --bfile cattle --cow --recode vcf --out cattle_vcf
r plink --vcf cattle_vcf.vcf --cow --make-bed --out from_vcf_fail
r plink --vcf cattle_vcf.vcf --cow --double-id --make-bed --out from_vcf
r plink --vcf cattle_vcf.vcf --cow --const-fid --make-bed --out from_vcf3
r plink --vcf rec_vcf-iid.vcf --cow --double-id --make-bed --out from_vcfiid
awk '$1=="NDA" {print $1,$2}' cattle.fam > nda_unrel.txt
r plink --bfile cattle --cow --keep nda_unrel.txt --make-bed --out nda_only
r plink --bfile cattle --cow --remove nda_unrel.txt --make-bed --out no_nda_unrel
r plink --bfile cattle --cow --chr 14 --from-mb 23 --to-mb 27 --make-bed --out chr14_plag1
r plink --bfile cattle --cow --chr 1-5,X --make-bed --out chr1to5X
r plink --bfile cattle --cow --snp snp02145 --window 500 --make-bed --out around_plag1
r plink --bfile cattle --cow --snps snp00010-snp00020 --make-bed --out snp_range
awk '$1==6 {print $2}' cattle.bim > chr6_snps.txt
r plink --bfile cattle --cow --extract chr6_snps.txt --make-bed --out chr6
r plink --bfile cattle --cow --exclude chr6_snps.txt --make-bed --out no_chr6
r plink --bfile cattle --cow --filter-females --make-bed --out females
r plink --bfile cattle --cow --allow-no-sex --filter-cases --make-bed --out cases_only
r plink --bfile cattle --cow --thin 0.25 --seed 7 --make-bed --out thin25
r plink --bfile cattle --cow --thin-count 1000 --seed 7 --make-bed --out thin1000
awk '$1 ~ /_F/ {print $1}' cattle.fam | sort -u > family_ids.txt
r plink --bfile cattle --cow --keep-fam family_ids.txt --make-bed --out families
r plink --bfile cattle --cow --autosome --make-bed --out autosomes
r plink --bfile cattle --cow --filter $D/breed.clst HOL --mfilter 1 --make-bed --out hol_filter
r plink --file $D/batch2 --cow --make-bed --out batch2
r plink --bfile cattle --cow --allow-no-sex --bmerge batch2 --make-bed --out merged_try1
r plink --bfile batch2 --cow --flip merged_try1-merge.missnp --make-bed --out batch2_flipped
r plink --bfile cattle --cow --bmerge batch2_flipped --make-bed --out merged_nosexflag
r plink --bfile cattle --cow --allow-no-sex --bmerge batch2_flipped --make-bed --out merged
awk '{print $1,$2,($1 ~ /B2$/)?2:1}' merged.fam > batch_as_pheno.txt
r plink --bfile merged --cow --allow-no-sex --pheno batch_as_pheno.txt --flip-scan --out merged_flipscan
python3 - <<'PY'
import random; random.seed(3)
lines=open('rec_ped.ped').read().split('\n')[:10]
out=[]
for l in lines:
    f=l.split(); g=f[6:]
    for j in range(0,len(g),2):
        if random.random()<0.01 and g[j]!='0': g[j],g[j+1]=g[j+1],g[j+1]
    out.append(' '.join(f[:6]+g))
open('regeno.ped','w').write('\n'.join(out)+'\n')
PY
cp rec_ped.map regeno.map
r plink --file regeno --cow --make-bed --out regeno
r plink --bfile cattle --cow --allow-no-sex --bmerge regeno --merge-mode 6 --out mm6
r plink --bfile cattle --cow --allow-no-sex --bmerge regeno --merge-mode 7 --out mm7
r plink --bfile cattle --cow --chr 1-10 --make-bed --out part1
r plink --bfile cattle --cow --chr 11-20 --make-bed --out part2
r plink --bfile cattle --cow --chr 21-29,X --make-bed --out part3
printf "part2\npart3\n" > mergelist.txt
r plink --bfile part1 --cow --allow-no-sex --merge-list mergelist.txt --make-bed --out rejoined
awk '{print $1,$2}' cattle.fam | sort -k2,2r > order_desc.txt
r plink --bfile cattle --cow --allow-no-sex --indiv-sort f order_desc.txt --make-bed --out resorted
echo "NDA NDA_022 1" > sexfix.txt
r plink --bfile cattle --cow --update-sex sexfix.txt --make-bed --out sexfixed
echo "HOL HOL_050_rep HOL HOL_050_B" > newids.txt
r plink --bfile cattle --cow --allow-no-sex --update-ids newids.txt --make-bed --out renamed
echo "snp02145 PLAG1_lead" > newname.txt
r plink --bfile cattle --cow --allow-no-sex --update-name newname.txt --make-bed --out snprenamed
awk '{print $2, $6}' cattle.bim > a2_as_ref.txt
r plink --bfile cattle --cow --allow-no-sex --a1-allele a2_as_ref.txt 2 1 --make-bed --out a1swapped
} > ch1_all.txt 2>&1
plink2 --bfile cattle --cow --make-pgen --out cattle_p2 > /dev/null 2>&1; plink2 --bfile cattle --cow --make-pgen vzs --out cattle_p2z > /dev/null 2>&1
plink2 --bfile cattle --cow --export vcf bgz --out cattle_p2vcf >/dev/null 2>&1
