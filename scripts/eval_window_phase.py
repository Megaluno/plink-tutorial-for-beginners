#!/usr/bin/env python3
"""Score PLINK 1.07 --hap-phase windows and a chromosome-wide phased VCF against the truth.
For every window of chip SNPs and every animal with at least two heterozygous SNPs in it,
the phase is correct if the two inferred haplotypes equal the two true ones.
Usage: eval_window_phase.py PREFIX(.phase-WIN*) MAP TRUTH.vcf.gz PHASED.vcf.gz"""
import sys, gzip, glob, re
import numpy as np
prefix, mapf, truthf, phasedf = sys.argv[1:5]
snps = [l.split()[1] for l in open(mapf)]
def read(fn):
    ids, d = None, {}
    for l in gzip.open(fn, "rt"):
        if l.startswith("##"): continue
        t = l.rstrip("\n").split("\t")
        if l.startswith("#"): ids = t[9:]; continue
        d[t[2]] = t
    return ids, d
tids, tru = read(truthf); pids, ph = read(phasedf)
def letters(t, k, which):
    a = t[9 + k].split(":")[0]; al = (t[3], t[4])
    return al[int(a[0 if which == 0 else 2])]
n_ok = n_all = m_ok = 0
files = glob.glob(prefix + ".phase-WIN*")
for fn in files:
    w = int(re.search(r"WIN(\d+)$", fn).group(1))
    win = snps[w - 1: w - 1 + 6]
    best = {}
    for line in open(fn):
        t = line.split()
        if t[0] == "FID" or len(t) < 7: continue
        if t[6] == "1": best[t[1]] = (t[3], t[4])
    for k, iid in enumerate(tids):
        h1 = "".join(letters(tru[s], k, 0) for s in win)
        h2 = "".join(letters(tru[s], k, 1) for s in win)
        nhet = sum(a != b for a, b in zip(h1, h2))
        if nhet < 2 or iid not in best: continue
        n_all += 1
        n_ok += {best[iid][0], best[iid][1]} == {h1, h2}
        kk = pids.index(iid)
        p1 = "".join(letters(ph[s], kk, 0) for s in win); p2 = "".join(letters(ph[s], kk, 1) for s in win)
        m_ok += {p1, p2} == {h1, h2}
print(f"{len(files)} windows; animal-windows with >=2 hets: {n_all}")
print(f"PLINK 1.07 --hap-phase correct: {n_ok / n_all:.3f}")
print(f"chromosome-wide Li-Stephens phase correct in the same windows: {m_ok / n_all:.3f}")
