#!/usr/bin/env python3
"""A small mixed-model GWAS (EMMAX-style, leave-one-chromosome-out), written for teaching.

Model: y = X b + g * beta + u + e,  u ~ N(0, s2g K),  e ~ N(0, s2e I)
K is the genomic relationship matrix (GRM) built from all chromosomes except the one being
tested (LOCO), so a SNP is never fitted twice. For each chromosome the variance ratio
delta = s2e / s2g is estimated once by REML (null model, no SNP), then every SNP is tested by
generalised least squares with delta fixed (the EMMAX / GCTA-mlma approximation).

Input: a PLINK .raw file from `plink --recode A` (counts of A1), the .bim file, a phenotype
file and a covariate file (FID IID columns, header).
Usage: lmm_gwas.py PREFIX.raw PREFIX.bim pheno.txt TRAIT covar.txt COV1,COV2,... OUT
"""
import sys
import numpy as np
from scipy import optimize, stats

raw, bim, phf, trait, covf, covnames, out = sys.argv[1:8]
covnames = [c for c in covnames.split(",") if c]

hdr = open(raw).readline().split()
dat = np.genfromtxt(raw, skip_header=1, dtype=str)
ids = [tuple(r[:2]) for r in dat]
G = np.where(dat[:, 6:] == "NA", "nan", dat[:, 6:]).astype(float)   # animals x SNPs, NA -> nan
snps = [h.rsplit("_", 1)[0] for h in hdr[6:]]
a1 = [h.rsplit("_", 1)[1] for h in hdr[6:]]
chrom = {l.split()[1]: l.split()[0] for l in open(bim)}
ch = np.array([chrom[s] for s in snps])

def table(fn, cols):
    h = open(fn).readline().split()
    d = {}
    for l in open(fn).readlines()[1:]:
        t = l.split(); d[(t[0], t[1])] = [float(t[h.index(c)]) for c in cols]
    return d
ph = table(phf, [trait]); cv = table(covf, covnames) if covnames else {}
keep = [k for k, i in enumerate(ids) if i in ph and ph[i][0] != -9 and (not covnames or i in cv)]
y = np.array([ph[ids[k]][0] for k in keep])
X = np.column_stack([np.ones(len(keep))] + ([np.array([cv[ids[k]] for k in keep])] if covnames else []))
G = G[keep]
n = len(y)

# standardised genotypes (mean-impute missing), as in the GRM of Chapter 10
p = np.nanmean(G, 0) / 2
ok = (p > 0) & (p < 1)
Z = np.where(np.isnan(G), 2 * p, G)
Zs = (Z - 2 * p) / np.sqrt(2 * p * (1 - p) + 1e-12)
auto = np.array([c.isdigit() for c in ch])

def reml_negll(logd, yr, Xr, S):
    d = np.exp(logd); w = 1 / (S + d)
    XtWX = Xr.T @ (Xr * w[:, None]); XtWy = Xr.T @ (yr * w)
    b = np.linalg.solve(XtWX, XtWy); r = yr - Xr @ b
    q = X.shape[1]
    s2 = (r * r * w).sum() / (n - q)
    return 0.5 * ((n - q) * np.log(s2) + np.log(S + d).sum() + np.linalg.slogdet(XtWX)[1])

res = []
h2s = []
for c in sorted(set(ch[auto]), key=int):
    other = ok & auto & (ch != c)
    K = Zs[:, other] @ Zs[:, other].T / other.sum()
    S, U = np.linalg.eigh(K)
    S = np.clip(S, 0, None)
    yr, Xr = U.T @ y, U.T @ X
    opt = optimize.minimize_scalar(reml_negll, bounds=(-6, 6), args=(yr, Xr, S), method="bounded")
    d = np.exp(opt.x); w = 1 / (S + d)
    h2s.append(1 / (1 + d))
    for j in np.where((ch == c) & ok)[0]:
        gr = U.T @ Z[:, j]
        Xa = np.column_stack([Xr, gr])
        XtWX = Xa.T @ (Xa * w[:, None]); XtWy = Xa.T @ (yr * w)
        b = np.linalg.solve(XtWX, XtWy); r = yr - Xa @ b
        s2 = (r * r * w).sum() / (n - Xa.shape[1])
        se = np.sqrt(s2 * np.linalg.inv(XtWX)[-1, -1])
        t = b[-1] / se
        res.append((c, snps[j], a1[j], b[-1], se, 2 * stats.t.sf(abs(t), n - Xa.shape[1])))

pv = np.array([r[5] for r in res])
lam = np.median(stats.chi2.isf(pv, 1)) / stats.chi2.ppf(0.5, 1)
with open(out, "w") as f:
    f.write("CHR\tSNP\tA1\tBETA\tSE\tP\n")
    for r in res: f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]:.4f}\t{r[5]:.4g}\n")
print(f"n = {n}, SNPs tested = {len(res)}, covariates = {covnames or 'none'}")
print(f"REML heritability (SNP-based, averaged over LOCO GRMs) = {np.mean(h2s):.2f}")
print(f"lambda = {lam:.3f}")
order = np.argsort(pv)
for k in order[:3]:
    r = res[k]; print(f"  {r[1]}  A1={r[2]}  BETA={r[3]:.2f}  P={r[5]:.3g}")
