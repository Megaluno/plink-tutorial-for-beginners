#!/usr/bin/env python3
"""Simulate data for the phasing and imputation chapter.

A coalescent simulation (msprime) of one 10 Mb stretch of cattle chromosome 14 in two
breeds that split 60 generations ago, each with an effective size of 150 since then.
It writes:
  ref_A_hd.vcf.gz      150 animals of breed A, high-density SNPs, PHASED (the reference panel)
  ref_B_hd.vcf.gz      150 animals of breed B, same SNPs, PHASED (a reference of another breed)
  target_chip.vcf.gz   100 other animals of breed A, chip SNPs only, UNPHASED, with some errors
  truth_target_hd.vcf.gz  the true phased haplotypes of the 100 targets at all HD SNPs
                       (answer key: never available in real life)
  target_pheno.txt     a weight trait with one QTL at a SNP that is NOT on the chip
  chr14.map            PLINK-format genetic map (1 cM per Mb) for the HD SNPs
Usage: python3 scripts/simulate_phasing.py [outdir]   (default data/phasing)
"""
import sys, os, gzip, json
import numpy as np
import msprime

OUT = sys.argv[1] if len(sys.argv) > 1 else "data/phasing"
os.makedirs(OUT, exist_ok=True)
SEED = 2026
rng = np.random.default_rng(SEED)
CHR, OFFSET, L = "14", 20_000_000, 10_000_000

dem = msprime.Demography()
dem.add_population(name="A", initial_size=150)
dem.add_population(name="B", initial_size=150)
dem.add_population(name="ANC", initial_size=4000)
dem.add_population_split(time=60, derived=["A", "B"], ancestral="ANC")
ts = msprime.sim_ancestry(samples={"A": 250, "B": 150}, demography=dem, sequence_length=L,
                          recombination_rate=1e-8, random_seed=SEED)
ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=SEED, model=msprime.BinaryMutationModel(),
                           discrete_genome=True)

G = ts.genotype_matrix()                     # sites x haplotypes (0/1), haplotypes paired per animal
pos = ts.tables.sites.position.astype(int)
popA = ts.samples(population=0); popB = ts.samples(population=1)
GA, GB = G[:, popA], G[:, popB]
refA, tgt = GA[:, :300], GA[:, 300:500]      # 150 reference animals, 100 target animals
fA = refA.mean(1)
keep = (np.minimum(fA, 1 - fA) >= 0.01) & (np.diff(np.r_[-1, pos]) > 0)
idx = np.where(keep)[0]
# thin to about one SNP per 3 kb: the "HD" panel
hd = []
last = -1e9
for i in idx:
    if pos[i] - last >= 3000:
        hd.append(i); last = pos[i]
hd = np.array(hd)
# chip: SNPs with MAF >= 0.05, about one per 30 kb, chosen from the HD panel
chip, last = [], -1e9
mafA = np.minimum(fA, 1 - fA)
for i in hd:
    if mafA[i] >= 0.05 and pos[i] - last >= 30000:
        chip.append(i); last = pos[i]
chip = np.array(chip)

BASES = "ACGT"
def alleles(n):
    out = []
    for _ in range(n):
        while True:
            a = rng.choice(4, 2, replace=False)
            if {BASES[a[0]], BASES[a[1]]} not in ({"A", "T"}, {"C", "G"}): break
        out.append((BASES[a[0]], BASES[a[1]]))
    return out
AL = dict(zip(hd, alleles(len(hd))))        # (REF, ALT); allele 0 = REF
ID = {i: f"hd{k + 1:05d}" for k, i in enumerate(hd)}

def write_vcf(fn, H, ids, sites, phased=True, gts=None):
    with gzip.open(fn, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CHR},length={OFFSET + L}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        for r, i in enumerate(sites):
            ref, alt = AL[i]
            if gts is not None:
                cells = gts[r]
            else:
                h = H[i]
                sep = "|" if phased else "/"
                cells = [f"{h[2*k]}{sep}{h[2*k+1]}" for k in range(len(ids))]
            f.write(f"{CHR}\t{OFFSET + pos[i]}\t{ID[i]}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(cells) + "\n")

idsA = [f"A_{k + 1:03d}" for k in range(150)]
idsB = [f"B_{k + 1:03d}" for k in range(150)]
idsT = [f"T_{k + 1:03d}" for k in range(100)]
write_vcf(f"{OUT}/ref_A_hd.vcf.gz", refA, idsA, hd)
write_vcf(f"{OUT}/ref_B_hd.vcf.gz", GB, idsB, hd)
write_vcf(f"{OUT}/truth_target_hd.vcf.gz", tgt, idsT, hd)

# target chip genotypes: unphased, allele order shuffled within animal, 0.5% missing, 0.1% errors
gts, nerr, nmiss = [], 0, 0
for i in chip:
    h = tgt[i]; row = []
    for k in range(100):
        g = int(h[2*k] + h[2*k+1])
        u = rng.random()
        if u < 0.005:
            row.append("./."); nmiss += 1; continue
        if u < 0.006:
            g = int(rng.choice([x for x in (0, 1, 2) if x != g])); nerr += 1
        row.append({0: "0/0", 1: "0/1", 2: "1/1"}[g])
    gts.append(row)
write_vcf(f"{OUT}/target_chip.vcf.gz", None, idsT, chip, gts=gts)

with open(f"{OUT}/chr14.map", "w") as f:       # PLINK map: chr, id, cM, bp (1 cM per Mb)
    for i in hd:
        bp = OFFSET + pos[i]
        f.write(f"{CHR}\t{ID[i]}\t{(bp - OFFSET) / 1e6:.6f}\t{bp}\n")

# phenotype: QTL at an HD SNP that is not on the chip, MAF ~0.3, +10 kg per ALT allele
not_chip = np.setdiff1d(hd, chip)
cand = [i for i in not_chip if 0.25 <= tgt[i].mean() <= 0.35 and 0.25 <= fA[i] <= 0.35]
qtl = cand[len(cand) // 2]
dos = tgt[qtl].reshape(-1, 2).sum(1)
y = 250 + 10 * dos + rng.normal(0, 10, 100)
with open(f"{OUT}/target_pheno.txt", "w") as f:
    f.write("IID\tWWT\n")
    for k in range(100):
        f.write(f"{idsT[k]}\t{y[k]:.1f}\n")

truth = dict(n_hd=int(len(hd)), n_chip=int(len(chip)), qtl=ID[qtl], qtl_bp=int(OFFSET + pos[qtl]),
             qtl_alt=AL[qtl][1], qtl_effect_per_alt="+10 kg", chip_missing=nmiss, chip_errors=nerr)
json.dump(truth, open(f"{OUT}/TRUTH_phasing.json", "w"), indent=1)
print(json.dumps(truth, indent=1))
