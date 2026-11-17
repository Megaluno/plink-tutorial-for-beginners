#!/usr/bin/env python3
"""Compare an imputed VCF (GT phased + DS) with the true haplotypes.
Usage: python3 eval_imputation.py IMPUTED.vcf.gz TRUTH.vcf.gz CHIP.vcf.gz [label]
Prints: genotype concordance and dosage r2 (imputed SNPs only) by MAF class,
the share of imputed SNPs with DR2 >= 0.8, and the switch error rate of the phase."""
import sys, gzip
import numpy as np

def read(fn, field="GT"):
    ids, data = None, {}
    with gzip.open(fn, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            t = line.rstrip("\n").split("\t")
            if line.startswith("#"): ids = t[9:]; continue
            data[t[2]] = t
    return ids, data

imp_ids, imp = read(sys.argv[1]); tru_ids, tru = read(sys.argv[2]); _, chip = read(sys.argv[3])
label = sys.argv[4] if len(sys.argv) > 4 else ""
n = len(imp_ids); col = [tru_ids.index(i) for i in imp_ids]
rows = []
for sid, t in imp.items():
    fmt = t[8].split(":"); ids = fmt.index("DS")
    cells = t[9:]
    ds = np.array([float(c.split(":")[ids]) for c in cells])
    gt = [c.split(":")[0] for c in cells]
    tc = [tru[sid][9 + k] for k in col]
    th = np.array([[int(c[0]), int(c[2])] for c in tc])
    g = th.sum(1)
    ih = np.array([[int(x[0]), int(x[2])] for x in gt])
    info = dict(kv.split("=") for kv in t[7].split(";") if "=" in kv)
    rows.append(dict(id=sid, typed=sid in chip, ds=ds, g=g, th=th, ih=ih,
                     dr2=float(info["DR2"]), pos=int(t[1])))

def maf(r):
    p = r["g"].mean() / 2
    return min(p, 1 - p)

imputed = [r for r in rows if not r["typed"]]
print(f"{label}: {n} animals, {len(imputed)} imputed SNPs")
print(f"{'MAF class':<12}{'SNPs':>6}{'concordance':>13}{'mean r2':>9}{'DR2>=0.8':>10}")
bins = [(0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.51)]
tot_c, tot_n = 0, 0
allr2, alldr2 = [], []
for lo, hi in bins + [(0, 0.51)]:
    sel = [r for r in imputed if lo <= maf(r) < hi]
    if not sel: continue
    conc = np.mean([np.mean(np.rint(r["ds"]) == r["g"]) for r in sel])
    r2 = []
    for r in sel:
        if r["g"].std() > 0 and r["ds"].std() > 0:
            r2.append(np.corrcoef(r["ds"], r["g"])[0, 1] ** 2)
        if (lo, hi) == (0, 0.51):
            if r["g"].std() > 0 and r["ds"].std() > 0:
                allr2.append(r2[-1]); alldr2.append(r["dr2"])
    good = np.mean([r["dr2"] >= 0.8 for r in sel])
    name = "all" if (lo, hi) == (0, 0.51) else f"{lo:.2f}-{min(hi, 0.5):.2f}"
    print(f"{name:<12}{len(sel):>6}{conc:>13.3f}{np.mean(r2):>9.3f}{good:>10.2f}")
print(f"correlation between DR2 and true r2 across SNPs: {np.corrcoef(allr2, alldr2)[0, 1]:.2f}")

# switch errors: at chip SNPs where the animal is truly heterozygous, compare phase orientation
typed = sorted([r for r in rows if r["typed"]], key=lambda r: r["pos"])
sw, opp = 0, 0
for k in range(n):
    prev = None
    for r in typed:
        a, b = r["th"][k]; x, y = r["ih"][k]
        if a == b or x == y: continue                 # need true het and inferred het
        same = (x == a)
        if prev is not None:
            opp += 1; sw += (same != prev)
        prev = same
print(f"switch error rate: {sw} switches in {opp} consecutive heterozygous pairs = {sw / opp:.3f}")
