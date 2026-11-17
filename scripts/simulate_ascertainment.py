#!/usr/bin/env python3
"""Simulate data for the ascertainment-bias section (Chapter 18).

Three breeds from one coalescent history (msprime):
  EUR1, EUR2  two European-type breeds that split 40 generations ago (Ne 150 each since),
              from a European ancestor with Ne 2,000, which split from AFR 1,500 generations ago;
  AFR         an African-type breed with a larger effective size (Ne 3,000), so more diversity.
Five 20 Mb chromosomes; 50 study animals per breed, plus 10 extra animals per breed used only
as "discovery panels" to choose chip SNPs, as array makers do.

Writes PLINK text files (.ped/.map) for four SNP sets typed on the same 150 study animals:
  seq        a random sample of all SNPs segregating in the 150 animals (stands in for sequence)
  chipEUR    SNPs with MAF >= 0.05 in a discovery panel of 10 EUR1 animals, about 1 per 25 kb
  chipALL    SNPs with MAF >= 0.05 in a discovery panel of 10 animals of EACH breed, same density
plus breeds.clst. Usage: python3 simulate_ascertainment.py [outdir]
"""
import sys, os
import numpy as np
import msprime

OUT = sys.argv[1] if len(sys.argv) > 1 else "data/ascertainment"
os.makedirs(OUT, exist_ok=True)
SEED = 2026
rng = np.random.default_rng(SEED)
NSTUDY, NDISC, NCHR, L = 50, 10, 5, 20_000_000
BREEDS = ["EUR1", "EUR2", "AFR"]

dem = msprime.Demography()
dem.add_population(name="EUR1", initial_size=150)
dem.add_population(name="EUR2", initial_size=150)
dem.add_population(name="AFR", initial_size=3000)
dem.add_population(name="EUR", initial_size=2000)
dem.add_population(name="ANC", initial_size=8000)
dem.add_population_split(time=40, derived=["EUR1", "EUR2"], ancestral="EUR")
dem.add_population_split(time=1500, derived=["EUR", "AFR"], ancestral="ANC")

sets = {k: [] for k in ("seq", "chipEUR", "chipALL")}   # list of (chrom, pos, genotypes 0/1/2 for study animals)
for c in range(1, NCHR + 1):
    ts = msprime.sim_ancestry(samples={b: NSTUDY + NDISC for b in BREEDS}, demography=dem, sequence_length=L,
                              recombination_rate=1e-8, random_seed=SEED + c)
    ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=SEED + c, model=msprime.BinaryMutationModel(),
                               discrete_genome=True)
    G = ts.genotype_matrix()
    pos = ts.tables.sites.position.astype(int)
    keep = np.r_[True, np.diff(pos) > 0]
    G, pos = G[keep], pos[keep]
    study, disc = [], {}
    for k, b in enumerate(BREEDS):
        hap = ts.samples(population=k)
        study.append(hap[: 2 * NSTUDY]); disc[b] = hap[2 * NSTUDY:]
    study = np.concatenate(study)
    Hs = G[:, study]
    geno = Hs[:, 0::2] + Hs[:, 1::2]                     # sites x 150 animals
    poly = (geno.sum(1) > 0) & (geno.sum(1) < 2 * geno.shape[1])

    def maf(cols):
        f = G[:, cols].mean(1); return np.minimum(f, 1 - f)
    def pick(ok, spacing):
        idx, last = [], -1e12
        for i in np.where(ok)[0]:
            if pos[i] - last >= spacing:
                idx.append(i); last = pos[i]
        return idx
    chipE = pick(poly & (maf(disc["EUR1"]) >= 0.05), 25_000)
    chipA = pick(poly & (maf(np.concatenate([disc[b] for b in BREEDS])) >= 0.05) &
                 np.any([maf(disc[b]) >= 0.05 for b in BREEDS], axis=0), 25_000)
    seqi = np.sort(rng.choice(np.where(poly)[0], size=len(chipE), replace=False))
    for name, idx in (("seq", seqi), ("chipEUR", chipE), ("chipALL", chipA)):
        sets[name] += [(c, int(pos[i]), geno[i]) for i in idx]

ids = [(b, f"{b}_{k + 1:03d}") for b in BREEDS for k in range(NSTUDY)]
with open(f"{OUT}/breeds.clst", "w") as f:
    for fid, iid in ids: f.write(f"{fid} {iid} {fid}\n")
for name, rows in sets.items():
    with open(f"{OUT}/{name}.map", "w") as f:
        for k, (c, p, _) in enumerate(rows): f.write(f"{c}\t{name}_{k + 1:05d}\t0\t{p}\n")
    M = np.array([g for _, _, g in rows])               # snps x animals
    code = {0: "A A", 1: "A G", 2: "G G"}
    with open(f"{OUT}/{name}.ped", "w") as f:
        for a, (fid, iid) in enumerate(ids):
            f.write(f"{fid} {iid} 0 0 0 -9 " + " ".join(code[int(x)] for x in M[:, a]) + "\n")
    print(name, len(rows), "SNPs")
