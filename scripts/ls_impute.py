#!/usr/bin/env python3
"""A small, transparent Li-Stephens phasing and imputation program, written for teaching.

It does the same job as BEAGLE, IMPUTE5 or Minimac in miniature:
  * each target animal is modelled as a mosaic of two reference haplotypes (a hidden Markov
    model whose hidden states are pairs of reference haplotypes);
  * the chip genotypes are the observations; switches between reference haplotypes are
    recombinations, more likely over larger genetic distances;
  * forward-backward probabilities give, at every chip SNP, how likely each reference
    haplotype is to be copied; these are interpolated to the untyped SNPs to give dosages;
  * the Viterbi path gives the most likely pair of haplotypes = the phased genotypes.
Real programs are far faster (they cluster haplotypes, use the PBWT, or select a few
"neighbour" haplotypes per animal) and also update the reference with the targets.

Usage:
  python3 ls_impute.py REF.vcf.gz TARGET.vcf.gz MAP OUTPREFIX [--ne 600] [--err 0.002] [--max-targets N]
Writes OUTPREFIX.vcf.gz (all reference SNPs: GT phased, DS dosage, GP; INFO: AF, DR2, IMP)
"""
import sys, gzip, argparse, time
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ref"); ap.add_argument("target"); ap.add_argument("map"); ap.add_argument("out")
ap.add_argument("--ne", type=float, default=600.0, help="effective size used in the switch rate")
ap.add_argument("--err", type=float, default=0.002, help="genotype error rate in the emissions")
ap.add_argument("--max-targets", type=int, default=0)
ap.add_argument("--dump-path", default="", help="write the Viterbi copying path of every target animal to this file")
args = ap.parse_args()

def read_vcf(fn):
    meta, ids, rows = None, [], []
    with gzip.open(fn, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            t = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                ids = t[9:]; continue
            rows.append(t)
    return ids, rows

rids, rrows = read_vcf(args.ref)
tids, trows = read_vcf(args.target)
cm = {}
for line in open(args.map):
    c, i, g, bp = line.split()
    cm[i] = float(g)

# reference haplotypes, sites x 2n
R = np.array([[int(a) for c in r[9:] for a in (c[0], c[2])] for r in rrows], dtype=np.int8)
rinfo = [(r[0], r[1], r[2], r[3], r[4]) for r in rrows]
rpos = {r[2]: k for k, r in enumerate(rrows)}
K = R.shape[1]
# target genotypes at the SNPs shared with the reference (matched by ID and alleles)
typed, Gt = [], []
for r in trows:
    k = rpos.get(r[2])
    if k is None or (r[3], r[4]) != (rinfo[k][3], rinfo[k][4]):
        continue          # not in the reference, or alleles do not match: skip
    typed.append(k)
    Gt.append([-1 if c[0] == "." else int(c[0]) + int(c[2]) for c in r[9:]])
typed = np.array(typed); Gt = np.array(Gt, dtype=np.int8)   # typed sites x targets
order = np.argsort(typed); typed = typed[order]; Gt = Gt[order]
nT = Gt.shape[1] if not args.max_targets else min(args.max_targets, Gt.shape[1])
gpos = np.array([cm[r[2]] for r in rinfo])                # cM of every reference site
m = len(typed)
d = np.diff(gpos[typed]) / 100.0                           # Morgans between typed sites
rho = 1 - np.exp(-4 * args.ne * d / K)                     # switch probability per haplotype copy

Rt = R[typed].astype(np.float64)                           # typed sites x K
print(f"reference: {K // 2} animals ({K} haplotypes), {R.shape[0]} SNPs; "
      f"target: {nT} animals, {m} typed SNPs shared with the reference", flush=True)

def emis(g, j):
    """K x K emission matrix at typed site j for observed genotype g"""
    if g < 0:
        return None
    h = Rt[j]
    pair = h[:, None] + h[None, :]
    e = args.err
    return np.where(pair == g, 1 - e, e / 2)

def transition(A, r):
    """apply the diploid Li-Stephens transition to a K x K probability matrix"""
    s = A.sum()
    row = A.sum(1, keepdims=True); col = A.sum(0, keepdims=True)
    return (1 - r) ** 2 * A + (1 - r) * (r / K) * (row + col) + (r / K) ** 2 * s

def impute_one(t):
    g = Gt[:, t]
    # forward
    F = np.empty((m, K, K))
    A = np.full((K, K), 1.0 / (K * K))
    for j in range(m):
        if j: A = transition(A, rho[j - 1])
        E = emis(g[j], j)
        if E is not None: A = A * E
        A /= A.sum(); F[j] = A
    # backward and posterior marginals per haplotype copy
    B = np.ones((K, K)) / (K * K)
    P1 = np.empty((m, K)); P2 = np.empty((m, K))
    Q2 = {}                                   # (typed marker, site) -> P(both copies carry ALT)
    for j in range(m - 1, -1, -1):
        post = F[j] * B; post /= post.sum()
        P1[j] = post.sum(1); P2[j] = post.sum(0)
        sites = near[j]
        if len(sites):
            Hs = R[sites].astype(np.float64)          # sites x K
            q2 = np.einsum("sk,kl,sl->s", Hs, post, Hs)
            for s_, v in zip(sites, q2): Q2[(j, s_)] = v
        if j:
            E = emis(g[j], j)
            X = B if E is None else B * E
            B = transition(X, rho[j - 1]); B /= B.sum()
    # Viterbi for a phased path (store the winning move per state)
    V = np.full((K, K), 1.0 / (K * K))
    E = emis(g[0], 0)
    if E is not None: V = V * E
    V /= V.max()
    moves = []
    for j in range(1, m):
        r = rho[j - 1]
        stay = (1 - r) ** 2 * V
        rmax = V.max(1); rarg = V.argmax(1)
        cmax = V.max(0); carg = V.argmax(0)
        vall = V.max(); varg = np.unravel_index(V.argmax(), V.shape)
        cand = np.stack([stay,
                         np.broadcast_to(((1 - r) * r / K) * rmax[:, None], (K, K)),
                         np.broadcast_to(((1 - r) * r / K) * cmax[None, :], (K, K)),
                         np.full((K, K), (r / K) ** 2 * vall)])
        best = cand.argmax(0).astype(np.int8)
        V = cand.max(0)
        E = emis(g[j], j)
        if E is not None: V = V * E
        V /= V.max()
        moves.append((best, rarg, carg, varg))
    k1, k2 = np.unravel_index(V.argmax(), V.shape)
    path = [(k1, k2)]
    for j in range(m - 1, 0, -1):
        best, rarg, carg, varg = moves[j - 1]
        b = best[k1, k2]
        if b == 1: k2 = rarg[k1]
        elif b == 2: k1 = carg[k2]
        elif b == 3: k1, k2 = varg
        path.append((k1, k2))
    path = path[::-1]
    return P1, P2, Q2, path

nsite = R.shape[0]
DS = np.zeros((nsite, nT)); H1 = np.zeros((nsite, nT), np.int8); H2 = np.zeros((nsite, nT), np.int8)
GP = np.zeros((nsite, nT, 3))
t0 = time.time()
# interpolation weights for every reference site between flanking typed sites (by cM)
left = np.searchsorted(typed, np.arange(nsite), side="right") - 1
near = [[] for _ in range(m)]
for s_ in range(nsite):
    l = left[s_]
    if l < 0: near[0].append(s_)
    elif l >= m - 1: near[m - 1].append(s_)
    else: near[l].append(s_); near[l + 1].append(s_)
near = [np.array(x, dtype=int) for x in near]
for t in range(nT):
    P1, P2, Q2, path = impute_one(t)
    if args.dump_path:
        with open(args.dump_path, "w" if t == 0 else "a") as f:
            if t == 0: f.write("IID\tID\tPOS\tCM\tHAP1_FROM\tHAP2_FROM\tGENO\n")
            for j in range(m):
                f.write(f"{tids[t]}\t{rinfo[typed[j]][2]}\t{rinfo[typed[j]][1]}\t{gpos[typed[j]]:.4f}\t"
                        f"{rids[path[j][0] // 2]}_{path[j][0] % 2 + 1}\t{rids[path[j][1] // 2]}_{path[j][1] % 2 + 1}\t{Gt[j, t]}\n")
    for s in range(nsite):
        l = left[s]
        if l < 0: a1, a2 = P1[0], P2[0]; pk = path[0]; q2 = Q2[(0, s)]
        elif l >= m - 1: a1, a2 = P1[m - 1], P2[m - 1]; pk = path[m - 1]; q2 = Q2[(m - 1, s)]
        else:
            g0, g1 = gpos[typed[l]], gpos[typed[l + 1]]
            w = 0.5 if g1 == g0 else (gpos[s] - g0) / (g1 - g0)
            a1 = (1 - w) * P1[l] + w * P1[l + 1]; a2 = (1 - w) * P2[l] + w * P2[l + 1]
            q2 = (1 - w) * Q2[(l, s)] + w * Q2[(l + 1, s)]
            pk = path[l] if w < 0.5 else path[l + 1]
        p1 = a1 @ R[s]; p2 = a2 @ R[s]          # P(ALT) on each haplotype copy
        DS[s, t] = p1 + p2
        # genotype probabilities from the joint posterior of the two copies
        GP[s, t] = np.clip([1 - p1 - p2 + q2, p1 + p2 - 2 * q2, q2], 0, 1)
        H1[s, t] = R[s, pk[0]]; H2[s, t] = R[s, pk[1]]
    if (t + 1) % 20 == 0:
        print(f"  {t + 1} animals done ({time.time() - t0:.0f} s)", flush=True)

# observed genotypes at typed SNPs override the model where they are not missing
istyped = np.zeros(nsite, bool); istyped[typed] = True
af = DS.mean(1) / 2
# DR2: estimated squared correlation between dosage and true genotype (as in BEAGLE)
var_ds = DS.var(1)
dr2 = np.where((af > 0) & (af < 1), var_ds / np.maximum(2 * af * (1 - af), 1e-12), 0.0)
dr2 = np.clip(dr2, 0, 1)
with gzip.open(args.out + ".vcf.gz", "wt") as f:
    f.write("##fileformat=VCFv4.2\n##source=ls_impute.py (teaching Li-Stephens imputation)\n")
    f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT allele frequency">\n')
    f.write('##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-squared: estimated squared correlation between dosage and true genotype">\n')
    f.write('##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker (not on the chip)">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
    f.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="ALT dosage">\n')
    f.write('##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(tids[:nT]) + "\n")
    for s in range(nsite):
        c, p, i, ref, alt = rinfo[s]
        info = f"AF={af[s]:.4f};DR2={dr2[s]:.3f}" + ("" if istyped[s] else ";IMP")
        cells = [f"{H1[s, t]}|{H2[s, t]}:{DS[s, t]:.3f}:{GP[s, t, 0]:.3f},{GP[s, t, 1]:.3f},{GP[s, t, 2]:.3f}"
                 for t in range(nT)]
        f.write(f"{c}\t{p}\t{i}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:DS:GP\t" + "\t".join(cells) + "\n")
print(f"wrote {args.out}.vcf.gz: {nsite} SNPs ({nsite - m} imputed), {nT} animals, "
      f"{time.time() - t0:.0f} s", flush=True)
