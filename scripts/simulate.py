#!/usr/bin/env python3
"""Simulate a small cattle SNP dataset for the PLINK workbook.
Three breeds (Angus ANG, Holstein HOL, N'Dama NDA), 29 autosomes + X,
families (sire, dam, 2 calves), traits, and deliberately planted QC problems.
Everything here is SIMULATED; gene-region positions are only loosely modelled on real ones."""
import numpy as np, os, json
rng = np.random.default_rng(2026)
OUT = "data"; os.makedirs(OUT, exist_ok=True)

# approximate bovine chromosome lengths (Mb), chr1..29, X
LEN = [158,137,121,121,121,119,113,114,105,104,107,91,84,85,85,81,75,66,64,72,71,61,52,63,43,51,46,46,51,139]
CHRS = [str(i) for i in range(1,30)] + ["X"]
NPER = 150; NX = 250
DENSE = {"5": (60.0, 64.0, 60), "6": (37.0, 40.0, 60), "14": (23.0, 27.0, 80)}  # chr: (start Mb, end Mb, n)
BASES = np.array(list("ACGT"))
COMP = {"A":"T","T":"A","C":"G","G":"C"}

snps = []  # dicts: chr, name, pos, a1, a2
for c, L in zip(CHRS, LEN):
    n = NX if c == "X" else NPER
    pos = list(rng.uniform(0.2, L - 0.2, n) * 1e6)
    if c in DENSE:
        s, e, k = DENSE[c]
        pos += list(rng.uniform(s, e, k) * 1e6)
    pos = sorted(set(int(p) for p in pos))
    for p in pos:
        while True:
            a = rng.choice(4, 2, replace=False)
            if {BASES[a[0]], BASES[a[1]]} not in ({"A","T"}, {"C","G"}): break
        snps.append(dict(chr=c, pos=p, a1=BASES[a[0]], a2=BASES[a[1]]))
# names: BovineHD-style made-up IDs
for i, s in enumerate(snps):
    s["name"] = f"BTA{s['chr']}_{s['pos']}" if False else f"snp{i+1:05d}"
M = len(snps)
chr_of = np.array([s["chr"] for s in snps]); pos_of = np.array([s["pos"] for s in snps])
# make ~3% of SNPs strand-ambiguous (A/T or C/G) for the strand chapter
amb = rng.choice(M, int(0.03 * M), replace=False)
for i in amb:
    if rng.random() < 0.5: snps[i]["a1"], snps[i]["a2"] = "A", "T"
    else: snps[i]["a1"], snps[i]["a2"] = "C", "G"

# ancestral and breed allele frequencies (Balding-Nichols); allele "1" = a2 coding of freq p (prob of a1)
p0 = np.clip(rng.beta(0.9, 0.9, M), 0.02, 0.98)
FST = {"ANG": 0.10, "HOL": 0.10, "NDA": 0.18}
pb = {}
for b, F in FST.items():
    a = p0 * (1 - F) / F; bb = (1 - p0) * (1 - F) / F
    pb[b] = np.clip(rng.beta(a, bb), 0.0, 1.0)
# a few rare variants (MAF < 1%) everywhere
rare = rng.choice(np.where(chr_of != "X")[0], 60, replace=False)
for b in pb: pb[b][rare] = rng.uniform(0.0, 0.008, rare.size)

# haplotypes with distance-decaying LD: a latent Gaussian AR(1) process along each chromosome,
# thresholded at the breed allele frequency (correlation exp(-d / LAMBDA) between neighbouring SNPs)
LAMBDA = 1.2e6
sign = np.where(rng.random(M) < 0.5, 1.0, -1.0)  # random LD phase per SNP
from statistics import NormalDist
_nd = NormalDist()
qcache = {b: np.array([_nd.inv_cdf(min(max(p, 1e-9), 1 - 1e-9)) for p in pb[b]]) for b in pb}
chr_idx = {c: np.where(chr_of == c)[0] for c in CHRS}
rho_prev = np.ones(M)
for c in CHRS:
    idx = chr_idx[c]
    rho_prev[idx[1:]] = np.exp(-np.diff(pos_of[idx]) / LAMBDA)
    rho_prev[idx[0]] = 0.0

def new_haps(b, n):
    z = np.empty((n, M))
    e = rng.standard_normal((n, M))
    for j in range(M):
        r = rho_prev[j]
        z[:, j] = e[:, j] if r == 0 else r * z[:, j - 1] + np.sqrt(1 - r * r) * e[:, j]
    zz = z * sign[None, :]
    return (zz < qcache[b][None, :]).astype(np.int8)  # 1 = carries a1, with probability p

_hapbuf = {}
def new_hap(b):
    if not _hapbuf.get(b):
        _hapbuf[b] = list(new_haps(b, 200))
    return _hapbuf[b].pop()

def recomb(h1, h2):
    """gamete from two haplotypes with ~1 crossover per 100 Mb"""
    g = np.empty(M, np.int8)
    for c in CHRS:
        idx = chr_idx[c]; L = pos_of[idx[-1]]
        nco = rng.poisson(L / 1e8); cps = np.sort(rng.uniform(0, L, nco))
        start = rng.integers(2); src = np.full(idx.size, start)
        for cp in cps: src[pos_of[idx] > cp] = 1 - src[pos_of[idx] > cp]
        g[idx] = np.where(src == 0, h1[idx], h2[idx])
    return g

animals = []  # dict: fid, iid, pat, mat, sex(1 m, 2 f), breed, h1, h2, true_sex
def add(fid, iid, pat, mat, sex, breed, h1, h2):
    animals.append(dict(fid=fid, iid=iid, pat=pat, mat=mat, sex=sex, true_sex=sex, breed=breed, h1=h1, h2=h2))
    return animals[-1]

xi = chr_idx["X"]
def fix_male_x(a):
    if a["true_sex"] == 1: a["h2"][xi] = a["h1"][xi]

NUNREL = {"ANG": 80, "HOL": 80, "NDA": 80}
for b, n in NUNREL.items():
    for i in range(1, n + 1):
        sex = 1 if rng.random() < 0.35 else 2
        a = add(b, f"{b}_{i:03d}", "0", "0", sex, b, new_hap(b), new_hap(b)); fix_male_x(a)

# families: 15 HOL + 15 NDA, sire + dam + 2 calves
for b in ["HOL", "NDA"]:
    for f in range(1, 16):
        fid = f"{b}_F{f:02d}"
        s = add(fid, f"{fid}_S", "0", "0", 1, b, new_hap(b), new_hap(b)); fix_male_x(s)
        d = add(fid, f"{fid}_D", "0", "0", 2, b, new_hap(b), new_hap(b))
        for c in (1, 2):
            sex = 1 if rng.random() < 0.5 else 2
            hm = recomb(d["h1"], d["h2"])
            hp = recomb(s["h1"], s["h2"])
            if sex == 1: hp[xi] = hm[xi]          # son: X from dam only (haploid)
            else: hp[xi] = s["h1"][xi]            # daughter: sire's X
            add(fid, f"{fid}_C{c}", s["iid"], d["iid"], sex, b, hm, hp)

byid = {a["iid"]: a for a in animals}
truth = {}
# --- planted problems in the GENOTYPES ---
# 1. misattributed paternity: HOL_F03_C2 true sire is HOL_011
c = byid["HOL_F03_C2"]; real = byid["HOL_011"]
hp = recomb(real["h1"], real["h2"])
if c["true_sex"] == 1: hp[xi] = c["h1"][xi]
else: hp[xi] = real["h1"][xi]
c["h2"] = hp; truth["wrong_sire"] = "HOL_F03_C2 (recorded sire HOL_F03_S, true sire HOL_011)"
# 2. inbred animal: NDA_040 haplotypes identical on chr 1-8
a = byid["NDA_040"]
for cc in [str(i) for i in range(1, 9)]: a["h2"][chr_idx[cc]] = a["h1"][chr_idx[cc]]
truth["inbred"] = "NDA_040 (homozygous across chr 1-8)"
# 3. crossbred outlier: HOL_077 has one NDA haplotype
byid["HOL_077"]["h2"] = new_hap("NDA"); fix_male_x(byid["HOL_077"]); truth["crossbred"] = "HOL_077 (half N'Dama)"
# 4. mislabelled breed: ANG_066 is really a Holstein
byid["ANG_066"]["h1"] = new_hap("HOL"); byid["ANG_066"]["h2"] = new_hap("HOL"); fix_male_x(byid["ANG_066"])
truth["mislabelled"] = "ANG_066 (labelled Angus, genetically Holstein)"

# genotype matrix: count of a1 alleles
G = np.array([a["h1"] + a["h2"] for a in animals], np.int8)
N = len(animals)
iids = [a["iid"] for a in animals]
# 5. duplicate sample HOL_050 genotyped twice (0.3% discordance)
dup = G[iids.index("HOL_050")].copy(); flip = rng.random(M) < 0.003
dup[flip] = rng.integers(0, 3, flip.sum())
src = byid["HOL_050"]
animals.append(dict(fid="HOL", iid="HOL_050_rep", pat="0", mat="0", sex=src["sex"], true_sex=src["true_sex"], breed="HOL", h1=None, h2=None))
G = np.vstack([G, dup]); iids.append("HOL_050_rep"); truth["duplicate"] = "HOL_050 and HOL_050_rep"
# 6. contaminated DNA: ANG_030 mixed with ANG_031 -> excess heterozygosity
i30, i31 = iids.index("ANG_030"), iids.index("ANG_031")
alle = np.stack([animals[i30]["h1"], animals[i30]["h2"], animals[i31]["h1"], animals[i31]["h2"]])
pick = np.array([rng.choice(4, 2, replace=False) for _ in range(M)])
mixed = alle[pick[:, 0], np.arange(M)] + alle[pick[:, 1], np.arange(M)]
# contamination drives calls toward het
tot = alle.sum(0)  # number of a1 copies among the 4 alleles
G[i30] = np.where(tot == 0, 0, np.where(tot == 4, 2, 1))  # any mixture of alleles is called heterozygous
truth["contaminated"] = "ANG_030 (DNA mixed with ANG_031)"
N = len(animals)

# genotype calling: -1 = missing
calls = G.astype(np.int16)
# background missingness 0.4 % and error 0.05 %
calls[rng.random(calls.shape) < 0.004] = -1
err = rng.random(calls.shape) < 0.0005
calls[err] = rng.integers(0, 3, err.sum())
# 7. high-missing SNPs (25)
auto = np.where(chr_of != "X")[0]
bad_snp = rng.choice(auto, 25, replace=False)
for j in bad_snp:
    r = rng.uniform(0.06, 0.30); calls[rng.random(N) < r, j] = -1
truth["high_missing_snps"] = [snps[j]["name"] for j in bad_snp]
# 8. high-missing animals (6)
bad_ind = ["ANG_012", "ANG_045", "HOL_023", "NDA_017", "NDA_061", "HOL_F07_C1"]
for nm, r in zip(bad_ind, [0.06, 0.25, 0.12, 0.09, 0.18, 0.15]):
    k = iids.index(nm); calls[k, rng.random(M) < r] = -1
truth["high_missing_animals"] = bad_ind
# 9. HWE-failing SNPs: 60% of heterozygotes miscalled as homozygote
hwe_bad = rng.choice(np.setdiff1d(auto, bad_snp), 12, replace=False)
for j in hwe_bad:
    hets = np.where(calls[:, j] == 1)[0]; sel = hets[rng.random(hets.size) < 0.6]
    calls[sel, j] = rng.choice([0, 2], sel.size)
truth["hwe_bad_snps"] = [snps[j]["name"] for j in hwe_bad]
# 10. Mendel-error SNPs (random 4% genotype error)
me_bad = rng.choice(np.setdiff1d(auto, np.concatenate([bad_snp, hwe_bad])), 4, replace=False)
for j in me_bad:
    e = rng.random(N) < 0.06; calls[e, j] = rng.integers(0, 3, e.sum())
truth["mendel_bad_snps"] = [snps[j]["name"] for j in me_bad]

# --- PHENOTYPES ---
breed = np.array([a["breed"] for a in animals])
sex = np.array([a["true_sex"] for a in animals])
env = (rng.random(N) < 0.5).astype(int)  # 1 = feedlot, 0 = pasture
herd = np.array([f"H{rng.integers(1, 7)}" for _ in range(N)])
age = rng.normal(205, 12, N).round(0)  # age at weaning (days)
def snp_in(c, lo, hi):
    w = np.where((chr_of == c) & (pos_of > lo * 1e6) & (pos_of < hi * 1e6))[0]
    # pick a common one (MAF > 0.2 overall)
    f = G[:, w].mean(0) / 2; w = w[(f > 0.25) & (f < 0.75)]
    return int(w[len(w) // 2])
qtl_plag1 = snp_in("14", 24.5, 25.5); qtl_lcorl = snp_in("6", 38.0, 39.0)
gxe_snp = snp_in("10", 40, 60); brd_snp = snp_in("5", 61.5, 62.5)
epi1 = snp_in("2", 50, 80); epi2 = snp_in("20", 30, 50)
g = G.astype(float)
beff = np.where(breed == "ANG", 20, np.where(breed == "HOL", 10, -40))
wwt = 200 + beff + 9 * g[:, qtl_plag1] + 5 * g[:, qtl_lcorl] + 12 * env * g[:, gxe_snp] \
      + np.where(sex == 1, 15, 0) + 0.3 * (age - 205) + rng.normal(0, 16, N)
bw = 30 + 0.12 * (wwt - 200) + rng.normal(0, 3, N)
ywt = 1.8 * wwt + rng.normal(0, 30, N)
liab = 0.9 * g[:, brd_snp] + 1.4 * ((g[:, epi1] >= 1) & (g[:, epi2] >= 1)) \
       + np.where(breed == "HOL", 0.8, np.where(breed == "NDA", -0.6, 0)) + 0.5 * env + rng.normal(0, 1, N)
brd = np.where(liab > np.quantile(liab, 0.60), 2, 1)  # 2 = case, 1 = control
# missing phenotypes for a few
miss_ph = rng.choice(N, 6, replace=False)
wwt[miss_ph[:3]] = np.nan
# 11. SNPs whose missingness depends on case status (differential failure)
diff = rng.choice(np.setdiff1d(auto, np.concatenate([bad_snp, hwe_bad, me_bad, [brd_snp, epi1, epi2, qtl_plag1, qtl_lcorl, gxe_snp]])), 8, replace=False)
for j in diff:
    r = np.where(brd == 2, 0.14, 0.01); calls[rng.random(N) < r, j] = -1
truth["diff_missing_snps"] = [snps[j]["name"] for j in diff]

truth["causal"] = dict(PLAG1_like_WWT=snps[qtl_plag1]["name"], LCORL_like_WWT=snps[qtl_lcorl]["name"],
                       GxE_WWT_feedlot=snps[gxe_snp]["name"], BRD=snps[brd_snp]["name"],
                       epistasis_pair=[snps[epi1]["name"], snps[epi2]["name"]])

# --- recorded (pedigree) sex with planted errors ---
rec_sex = sex.copy()
k = iids.index("NDA_022"); rec_sex[k] = 2 if sex[k] == 1 else 1
k = iids.index("ANG_071"); rec_sex[k] = 2 if sex[k] == 1 else 1
k = iids.index("HOL_064"); rec_sex[k] = 0   # unknown sex
k = iids.index("NDA_009"); rec_sex[k] = 0
truth["sex_errors"] = {"NDA_022": "recorded sex flipped", "ANG_071": "recorded sex flipped", "HOL_064": "sex unknown (0)", "NDA_009": "sex unknown (0)"}
truth["true_sex"] = {nm: int(sex[iids.index(nm)]) for nm in ["NDA_022", "ANG_071", "HOL_064", "NDA_009"]}

# --- write PED/MAP (ACGT alleles) ---
a1 = np.array([s["a1"] for s in snps]); a2 = np.array([s["a2"] for s in snps])
def gt_strings(row):
    out = []
    for j, v in enumerate(row):
        if v < 0: out.append("0 0")
        elif v == 2: out.append(f"{a1[j]} {a1[j]}")
        elif v == 1: out.append(f"{a1[j]} {a2[j]}")
        else: out.append(f"{a2[j]} {a2[j]}")
    return " ".join(out)

with open(f"{OUT}/cattle.map", "w") as f:
    for s in snps:
        cm = round(s["pos"] / 1e6, 4)
        f.write(f"{s['chr']}\t{s['name']}\t{cm}\t{s['pos']}\n")
with open(f"{OUT}/cattle.ped", "w") as f:
    for k, a in enumerate(animals):
        f.write(f"{a['fid']} {a['iid']} {a['pat']} {a['mat']} {rec_sex[k]} {brd[k]} {gt_strings(calls[k])}\n")

def fmt(x): return "-9" if np.isnan(x) else f"{x:.1f}"
with open(f"{OUT}/pheno.txt", "w") as f:
    f.write("FID IID WWT BW YWT BRD\n")
    for k, a in enumerate(animals):
        f.write(f"{a['fid']} {a['iid']} {fmt(wwt[k])} {fmt(bw[k])} {fmt(ywt[k])} {brd[k]}\n")
with open(f"{OUT}/covar.txt", "w") as f:
    f.write("FID IID SEX_M FEEDLOT AGE HERD_N\n")
    for k, a in enumerate(animals):
        f.write(f"{a['fid']} {a['iid']} {int(sex[k]==1)} {env[k]} {int(age[k])} {herd[k][1:]}\n")
with open(f"{OUT}/env.txt", "w") as f:   # for --gxe: 1/2 coded groups
    for k, a in enumerate(animals): f.write(f"{a['fid']} {a['iid']} {env[k]+1}\n")
with open(f"{OUT}/breed.clst", "w") as f:  # cluster file for --within
    for k, a in enumerate(animals): f.write(f"{a['fid']} {a['iid']} {a['breed']}\n")

# --- second genotyping batch for merging (40 new animals, different SNP subset, some on opposite strand) ---
B2 = 40
b2b = np.repeat(["ANG", "HOL", "NDA", "NDA"], 10)
H = np.array([new_hap(b) + new_hap(b) for b in b2b], np.int8)
b2sex = np.where(rng.random(B2) < 0.4, 1, 2)
for k in range(B2):
    if b2sex[k] == 1: H[k, xi] = 2 * (H[k, xi] > 0)  # crude haploid X
keep2 = np.sort(rng.choice(M, int(0.85 * M), replace=False))  # batch 2 typed 85% of SNPs
extra = 120  # plus 120 SNPs not in batch 1
flipset = rng.choice(keep2, 40, replace=False)  # reported on the opposite strand
truth["batch2_flipped"] = [snps[j]["name"] for j in sorted(flipset)]
with open(f"{OUT}/batch2.map", "w") as fm:
    for j in keep2: fm.write(f"{snps[j]['chr']}\t{snps[j]['name']}\t0\t{snps[j]['pos']}\n")
    ex_pos = []
    for e in range(extra):
        c = rng.choice(CHRS[:29]); p = int(rng.uniform(1e6, LEN[int(c)-1] * 1e6)); ex_pos.append((c, p))
        fm.write(f"{c}\tb2snp{e+1:04d}\t0\t{p}\n")
ex_alle = [rng.choice(4, 2, replace=False) for _ in range(extra)]
with open(f"{OUT}/batch2.ped", "w") as fp:
    for k in range(B2):
        row = []
        for j in keep2:
            x, y = a1[j], a2[j]
            if j in flipset: x, y = COMP[x], COMP[y]
            v = H[k, j]
            if rng.random() < 0.004: row.append("0 0"); continue
            row.append(f"{x} {x}" if v == 2 else f"{x} {y}" if v == 1 else f"{y} {y}")
        for e in range(extra):
            u, w = BASES[ex_alle[e][0]], BASES[ex_alle[e][1]]; v = rng.binomial(2, 0.4)
            row.append(f"{u} {u}" if v == 2 else f"{u} {w}" if v == 1 else f"{w} {w}")
        fp.write(f"{b2b[k]}B2 {b2b[k]}_B2_{k+1:03d} 0 0 {b2sex[k]} -9 " + " ".join(row) + "\n")

# --- imputed dosage file for the chr14 dense region (fractional a1 counts) ---
w14 = np.where((chr_of == "14") & (pos_of > 23e6) & (pos_of < 27e6))[0]
with open(f"{OUT}/chr14.dosage", "w") as f:
    f.write("SNP A1 A2 " + " ".join(f"{a['fid']} {a['iid']}" for a in animals) + "\n")
    for j in w14:
        d = np.clip(G[:, j] + rng.normal(0, 0.15, N), 0, 2)
        f.write(f"{snps[j]['name']} {a1[j]} {a2[j]} " + " ".join(f"{x:.3f}" for x in d) + "\n")
with open(f"{OUT}/chr14.dosage.map", "w") as f:
    for j in w14: f.write(f"14\t{snps[j]['name']}\t0\t{snps[j]['pos']}\n")
with open(f"{OUT}/dosage.list", "w") as f: f.write("chr14.dosage\n")

# --- CNV calls (segmental, rare): .cnv file ---
cases = set(iids[k] for k in range(N) if brd[k] == 2)
cnv_hot = ("7", 40.0e6, 40.4e6)
with open(f"{OUT}/cattle.cnv", "w") as f:
    f.write("FID IID CHR BP1 BP2 TYPE SCORE SITES\n")
    for k, a in enumerate(animals):
        lam = 4.0 if a["iid"] in cases else 2.8
        for _ in range(rng.poisson(lam)):
            c = rng.integers(1, 30); L = LEN[c - 1] * 1e6
            s = int(rng.uniform(1e6, L - 2e6)); ln = int(rng.lognormal(np.log(80e3), 0.8))
            typ = 1 if rng.random() < 0.6 else 3
            f.write(f"{a['fid']} {a['iid']} {c} {s} {s+ln} {typ} {rng.uniform(5,40):.2f} {max(3, ln//8000)}\n")
        if (a["iid"] in cases and rng.random() < 0.22) or (a["iid"] not in cases and rng.random() < 0.04):
            s = int(cnv_hot[1] + rng.uniform(-20e3, 20e3)); e = int(cnv_hot[2] + rng.uniform(-20e3, 20e3))
            f.write(f"{a['fid']} {a['iid']} 7 {s} {e} 1 {rng.uniform(10,40):.2f} {(e-s)//8000}\n")
truth["cnv_hotspot"] = "chr7:40.0-40.4 Mb deletion enriched in BRD cases; cases carry ~40% more CNVs overall"

# gene / set file (illustrative regions, positions only loosely modelled on real genes)
genes = [("14", 24.0e6, 26.0e6, "PLAG1_region"), ("6", 37.5e6, 39.5e6, "LCORL_region"),
         ("5", 61.0e6, 63.0e6, "BRD_region"), ("10", 40e6, 60e6, "GXE_region"),
         ("1", 10e6, 12e6, "NULL_A"), ("3", 50e6, 52e6, "NULL_B"), ("8", 70e6, 72e6, "NULL_C"), ("12", 30e6, 32e6, "NULL_D")]
with open(f"{OUT}/genes.txt", "w") as f:
    for c, s, e, n in genes: f.write(f"{c} {int(s)} {int(e)} {n}\n")
json.dump(truth, open(f"{OUT}/TRUTH.json", "w"), indent=1)
print("animals", N, "snps", M, "batch2 snps", len(keep2) + extra)
print(json.dumps(truth["causal"]))
