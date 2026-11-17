"""Make every figure for the PLINK tutorial website from real PLINK output files."""
import os, json, collections
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

W = "work"      # folder that holds the PLINK outputs; change to your own
OUT = "images"  # where the figures are written
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito colour-blind safe palette
BLUE, ORANGE, GREEN, RED, PURPLE, SKY, YELLOW, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#7f7f7f"
BREED_COL = {"ANG": BLUE, "HOL": ORANGE, "NDA": GREEN}
BREED_NAME = {"ANG": "Angus", "HOL": "Holstein", "NDA": "N'Dama"}
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlesize": 12, "figure.dpi": 110,
                     "savefig.dpi": 150, "savefig.bbox": "tight"})

def save(name):
    plt.savefig(f"{OUT}/{name}.png")
    plt.close()

def ws(path, **kw):
    return pd.read_csv(path, sep=r"\s+", **kw)

breed = {l.split()[1]: l.split()[2] for l in open(f"{W}/data/breed.clst")}

# ---------------------------------------------------------------- 1. file formats schematic
def fig_formats():
    fig, ax = plt.subplots(figsize=(10, 5.2)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 60)
    def box(x, y, w, h, title, lines, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=col, ec="#333", alpha=0.18))
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc="none", ec="#333"))
        ax.text(x + 1.2, y + h - 2.6, title, fontweight="bold", fontsize=11)
        for i, t in enumerate(lines):
            ax.text(x + 1.2, y + h - 6.2 - i * 3.1, t, family="monospace", fontsize=8.6)
    box(1, 31, 44, 27, "cattle.ped  (one line per animal)",
        ["FID  IID      sire dam sex phen  genotypes",
         "ANG  ANG_001  0    0   2   1     A A  T C ...",
         "HOL_F01 HOL_F01_C1 HOL_F01_S HOL_F01_D 1 2 ...",
         "", "6 columns about the animal,", "then 2 letters per SNP (0 0 = missing)"], ORANGE)
    box(1, 3, 44, 22, "cattle.map  (one line per SNP)",
        ["chr  SNP       cM      bp", "1    snp00001  0.484   484050", "1    snp00002  2.2237  2223660", "",
         "same order as the genotype columns"], ORANGE)
    ax.annotate("", xy=(55, 30), xytext=(47, 30), arrowprops=dict(arrowstyle="-|>", lw=2, color="#333"))
    ax.text(51, 33.5, "--make-bed", ha="center", fontsize=10, family="monospace")
    box(56, 42, 43, 16, "cattle.fam  (= first 6 columns of .ped)", ["ANG ANG_001 0 0 2 1", "ANG ANG_002 0 0 1 2"], BLUE)
    box(56, 22, 43, 16, "cattle.bim  (= .map + the two alleles)", ["1 snp00001 0.484 484050 G A", "1 snp00002 2.2237 2223660 C A"], BLUE)
    box(56, 3, 43, 15, "cattle.bed  (genotypes, 2 bits each)", ["binary: do not open in a text editor", "6.9 MB .ped  ->  0.44 MB .bed"], BLUE)
    ax.text(1, 59.5, "Text format (easy to read, large)", color="#8a5a00", fontweight="bold")
    ax.text(56, 59.5, "Binary format (compact, fast)", color="#004f7c", fontweight="bold")
    save("01_file_formats")

# ---------------------------------------------------------------- 2. toy missingness grid
def toy_grid(ax, data, title, crossed_rows=(), crossed_cols=()):
    n, m = data.shape
    ax.set_xlim(-1.6, m + 2.2); ax.set_ylim(-1.4, n + 0.9); ax.axis("off"); ax.invert_yaxis()
    for i in range(n):
        for j in range(m):
            missing = data[i, j] in ("00", "--")
            fc = "#fde0d9" if missing else "#e8f1f8"
            if i in crossed_rows or j in crossed_cols: fc = "#eeeeee"
            ax.add_patch(Rectangle((j, i), 1, 1, fc=fc, ec="white", lw=2))
            ax.text(j + 0.5, i + 0.5, data[i, j], ha="center", va="center", family="monospace",
                    color=RED if missing else "#333", fontweight="bold" if missing else "normal", fontsize=10)
    for j in range(m):
        ax.text(j + 0.5, -0.45, f"SNP{j+1}", ha="center", fontsize=8.5)
        col = [data[i, j] for i in range(n) if i not in crossed_rows]
        rate = sum(c in ("00", "--") for c in col) / len(col)
        ax.text(j + 0.5, n + 0.45, f"{rate:.0%}", ha="center", fontsize=9, color=RED if rate > 0.1 else "#333")
    for i in range(n):
        ax.text(-0.25, i + 0.5, f"Animal {i+1}", ha="right", va="center", fontsize=8.5)
        row = [data[i, j] for j in range(m) if j not in crossed_cols]
        rate = sum(c in ("00", "--") for c in row) / len(row)
        ax.text(m + 0.2, i + 0.5, f"{rate:.0%}", va="center", fontsize=9, color=RED if rate > 0.1 else "#333")
    for i in crossed_rows: ax.plot([0, m], [i + 0.5, i + 0.5], color=RED, lw=2)
    for j in crossed_cols: ax.plot([j + 0.5, j + 0.5], [0, n], color=RED, lw=2)
    ax.text(m + 0.2, -0.45, "per animal", fontsize=8.5, style="italic")
    ax.text(-1.5, n + 0.45, "per SNP", fontsize=8.5, style="italic")
    ax.set_title(title, fontsize=11)

def fig_toy_missing():
    # Illumina-style A/B genotypes; "--" = missing
    # case 1: SNP3 failed in four of six animals
    a = np.array([["AA", "AB", "--", "BB", "AA", "AB", "AB", "BB"],
                  ["AB", "AA", "--", "AB", "BB", "AA", "--", "AB"],
                  ["BB", "AB", "--", "AA", "AB", "BB", "AA", "AB"],
                  ["AA", "BB", "--", "AB", "AA", "AB", "BB", "AA"],
                  ["AB", "AB", "AA", "BB", "AB", "AA", "AB", "BB"],
                  ["BB", "AA", "AB", "AB", "AA", "BB", "AB", "AA"]])
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
    toy_grid(axs[0], a, "Before: four of six animals look poor (12 to 25% missing)")
    toy_grid(axs[1], a, "After dropping SNP3 first (--geno 0.5)", crossed_cols=(2,))
    save("02_toy_bad_snp")
    # case 2: Animal 4's sample failed at six of eight SNPs
    b = np.array([["AA", "AB", "BB", "BB", "AA", "AB", "AB", "BB"],
                  ["AB", "AA", "AB", "AB", "BB", "AA", "AA", "AB"],
                  ["BB", "AB", "AA", "AA", "AB", "BB", "AA", "AB"],
                  ["--", "--", "AB", "--", "--", "--", "BB", "--"],
                  ["AB", "AB", "AA", "BB", "AB", "AA", "AB", "BB"],
                  ["BB", "AA", "AB", "AB", "AA", "BB", "AB", "AA"]])
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
    toy_grid(axs[0], b, "Before: six of eight SNPs look poor (17% missing)")
    toy_grid(axs[1], b, "After dropping Animal 4 first (--mind 0.5)", crossed_rows=(3,))
    save("02_toy_bad_animal")

# ---------------------------------------------------------------- 3. real missingness histograms
def fig_missing():
    im = ws(f"{W}/figdata/miss_raw.imiss"); lm = ws(f"{W}/figdata/miss_raw.lmiss")
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    axs[0].hist(lm.F_MISS * 100, bins=np.arange(0, 32, 0.5), color=BLUE)
    axs[0].set_yscale("log"); axs[0].axvline(5, color=RED, ls="--"); axs[0].text(5.5, 1500, "--geno 0.05", color=RED)
    axs[0].set_xlabel("Missing genotypes per SNP (%)"); axs[0].set_ylabel("Number of SNPs (log scale)")
    axs[0].set_title("Per SNP (.lmiss): 4,800 SNPs")
    axs[1].hist(im.F_MISS * 100, bins=np.arange(0, 26, 0.5), color=ORANGE)
    axs[1].set_yscale("log"); axs[1].axvline(5, color=RED, ls="--"); axs[1].text(5.5, 150, "--mind 0.05", color=RED)
    for _, r in im[im.F_MISS > 0.05].iterrows():
        axs[1].annotate(r.IID, (r.F_MISS * 100, 1), xytext=(r.F_MISS * 100, 3.5 + (hash(r.IID) % 4) * 3), fontsize=7.5,
                        ha="center", arrowprops=dict(arrowstyle="-", color=GREY, lw=0.6))
    axs[1].set_xlabel("Missing genotypes per animal (%)"); axs[1].set_ylabel("Number of animals (log scale)")
    axs[1].set_title("Per animal (.imiss): 361 animals")
    save("03_missingness_hist")

    het = ws(f"{W}/figdata/het_raw.het")
    d = het.merge(im, on=["FID", "IID"])
    d["HET"] = (d["N(NM)"] - d["O(HOM)"]) / d["N(NM)"]
    m, s = d.HET.mean(), d.HET.std()
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.scatter(d.F_MISS * 100, d.HET, s=14, color=GREY, alpha=0.6)
    ax.axhspan(m - 3 * s, m + 3 * s, color=GREEN, alpha=0.08); ax.axvline(5, color=RED, ls="--")
    ax.text(5.4, m + 3 * s + 0.004, "--mind 0.05", color=RED, fontsize=9)
    ax.text(16, m + 3 * s - 0.012, "mean ± 3 SD", color=GREEN, fontsize=9)
    for _, r in d[(d.F_MISS > 0.05) | (abs(d.HET - m) > 3 * s)].iterrows():
        ax.scatter(r.F_MISS * 100, r.HET, s=30, color=RED)
        ax.annotate(r.IID, (r.F_MISS * 100, r.HET), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Missing genotypes per animal (%)"); ax.set_ylabel("Observed heterozygosity")
    ax.set_title("Call rate against heterozygosity: the two sample checks together")
    save("03_miss_vs_het")

# ---------------------------------------------------------------- 4. MAF and HWE
def fig_maf_hwe():
    f = ws(f"{W}/ch2/freq.frq")
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.hist(f.MAF.dropna(), bins=np.arange(0, 0.51, 0.01), color=BLUE)
    for t, c in [(0.01, RED), (0.05, ORANGE)]:
        ax.axvline(t, color=c, ls="--"); ax.text(t + 0.004, ax.get_ylim()[1] * 0.85 if t == 0.01 else ax.get_ylim()[1] * 0.7, f"--maf {t}", color=c, fontsize=9)
    ax.set_xlabel("Minor allele frequency"); ax.set_ylabel("Number of SNPs"); ax.set_title("MAF of 4,800 SNPs (founders only)")
    save("04_maf_hist")

    truth = set(json.load(open(f"{W}/data/TRUTH.json"))["hwe_bad_snps"])
    h = ws(f"{W}/ch2/hwe.hwe"); h = h[h.TEST == "ALL"].dropna(subset=["P"])
    within = []
    for b in ["ANG", "HOL", "NDA"]:
        x = ws(f"{W}/ch2/hwe_{b}.hwe"); within.append(x[x.TEST == "ALL"][["SNP", "P"]].set_index("SNP").P)
    wmin = pd.concat(within, axis=1).min(axis=1)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, p, title in [(axs[0], h.set_index("SNP").P, "All breeds pooled"), (axs[1], wmin, "Within breed (smallest of 3 p-values)")]:
        p = p.dropna(); lp = -np.log10(p.clip(lower=1e-300)); bad = p.index.isin(truth)
        ax.scatter(np.arange(len(p))[~bad], lp[~bad], s=5, color=GREY, alpha=0.6, label="normal SNP")
        ax.scatter(np.arange(len(p))[bad], lp[bad], s=26, color=RED, label="planted calling error")
        ax.axhline(6, color=BLUE, ls="--"); ax.text(10, 6.5, "p = 1e-6", color=BLUE, fontsize=9)
        n = (p < 1e-6).sum(); ax.set_title(f"{title}: {n} SNPs below 1e-6")
        ax.set_xlabel("SNPs in genome order"); ax.set_ylabel("-log10(HWE p-value)"); ax.set_ylim(0, 40)
    axs[0].legend(loc="upper right", fontsize=8.5)
    save("04_hwe_pooled_vs_within")

    # worked HWE example bar chart
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    exp = [9, 42, 49]; obs = [20, 20, 60]; lab = ["AA", "AG", "GG"]; x = np.arange(3)
    ax.bar(x - 0.18, exp, 0.36, color=GREY, label="expected (p², 2pq, q²)")
    ax.bar(x + 0.18, obs, 0.36, color=RED, label="observed (too few heterozygotes)")
    for i in range(3):
        ax.text(x[i] - 0.18, exp[i] + 1, exp[i], ha="center", fontsize=9); ax.text(x[i] + 0.18, obs[i] + 1, obs[i], ha="center", fontsize=9)
    ax.set_xticks(x, lab); ax.set_ylabel("Animals (out of 100)"); ax.legend(fontsize=8.5, frameon=False)
    ax.set_title("Hardy-Weinberg: p(A) = 0.3, q(G) = 0.7")
    save("04_hwe_example")

# ---------------------------------------------------------------- 5. sex check and het F
def fig_sex_het():
    s = ws(f"{W}/ch2/sexcheck.sexcheck")
    fig, ax = plt.subplots(figsize=(7.5, 4))
    bins = np.arange(-0.3, 1.05, 0.025)
    for code, col, lab in [(1, BLUE, "recorded male"), (2, ORANGE, "recorded female"), (0, GREY, "sex unknown")]:
        ax.hist(s[s.PEDSEX == code].F, bins=bins, color=col, alpha=0.75, label=lab)
    for t in (0.2, 0.8): ax.axvline(t, color=RED, ls="--")
    ax.text(0.03, ax.get_ylim()[1] * 0.92, "called female", fontsize=9); ax.text(0.82, ax.get_ylim()[1] * 0.92, "called male", fontsize=9)
    ax.text(0.33, ax.get_ylim()[1] * 0.92, "'problem' zone", fontsize=9, color=RED)
    for _, r in s[(s.PEDSEX != s.SNPSEX) & (s.F > 0.8)].iterrows():
        ax.annotate(f"{r.IID} (recorded {r.PEDSEX})", (r.F, 1), xytext=(0.62, 25 + 8 * list(s.IID).index(r.IID) % 30), fontsize=7.5,
                    arrowprops=dict(arrowstyle="-", color=GREY, lw=0.6))
    ax.set_xlabel("F on the X chromosome"); ax.set_ylabel("Animals"); ax.legend(fontsize=8.5, loc="center left", frameon=False)
    ax.set_title("--check-sex: bulls sit near F = 1, cows near F = 0")
    save("05_checksex")

    h = ws(f"{W}/ch2/het.het"); h["breed"] = h.IID.map(breed)
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for b in ["ANG", "HOL", "NDA"]:
        ax.hist(h[h.breed == b].F, bins=np.arange(-0.8, 0.5, 0.02), color=BREED_COL[b], alpha=0.7, label=BREED_NAME[b])
    m, sd = h.F.mean(), h.F.std()
    for k in (-3, 3): ax.axvline(m + k * sd, color=RED, ls="--")
    for iid, txt in [("ANG_030", "ANG_030: contaminated\n(too heterozygous)"), ("NDA_040", "NDA_040: inbred")]:
        v = h[h.IID == iid].F.values[0]; ax.annotate(txt, (v, 0.5), xytext=(v + (0.05 if v < 0 else -0.25), 18), fontsize=8.5,
                                                     arrowprops=dict(arrowstyle="->", color="#333"))
    ax.set_xlabel("Inbreeding coefficient F (--het)"); ax.set_ylabel("Animals"); ax.legend(fontsize=8.5, frameon=False)
    ax.set_title("Heterozygosity check: red lines are mean ± 3 SD")
    save("05_het_F")

# ---------------------------------------------------------------- 6. IBD
def fig_ibd():
    g = ws(f"{W}/ch2/ibd.genome")
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    cols = {"UN": GREY, "OT": SKY, "FS": GREEN, "PO": ORANGE}
    labs = {"UN": "different family (UN)", "OT": "same FID, not parent/sib (OT)", "FS": "recorded full sibs (FS)", "PO": "recorded parent-offspring (PO)"}
    for rt in ["UN", "OT", "FS", "PO"]:
        x = g[g.RT == rt]
        if rt == "UN": x = x.sample(6000, random_state=1)
        ax.scatter(x.Z0, x.Z1, s=8 if rt in ("UN", "OT") else 22, color=cols[rt], alpha=0.5 if rt in ("UN", "OT") else 0.9, label=labs[rt])
    for (z0, z1, t) in [(0, 1, "parent-offspring\n(Z0=0, Z1=1)"), (0.25, 0.5, "full sibs\n(0.25, 0.5)"), (0.5, 0.5, "half sibs\n(0.5, 0.5)"), (1, 0, "unrelated\n(1, 0)"), (0, 0, "duplicate\n(0, 0)")]:
        ax.scatter(z0, z1, marker="X", s=90, color="black"); ax.annotate(t, (z0, z1), xytext=(6, 6), textcoords="offset points", fontsize=8)
    for a, b, t in [("HOL_F03_S", "HOL_F03_C2", "recorded sire\nbut Z0 = 1"), ("HOL_011", "HOL_F03_C2", "HOL_011: the\nreal sire")]:
        r = g[((g.IID1 == a) & (g.IID2 == b)) | ((g.IID1 == b) & (g.IID2 == a))].iloc[0]
        ax.scatter(r.Z0, r.Z1, s=70, facecolor="none", edgecolor=RED, lw=2)
        ax.annotate(t, (r.Z0, r.Z1), xytext=(-95, 45) if r.Z0 > 0.5 else (45, -55), textcoords="offset points", fontsize=8, color=RED,
                    arrowprops=dict(arrowstyle="->", color=RED))
    ax.set_xlabel("Z0 = P(share 0 alleles IBD)"); ax.set_ylabel("Z1 = P(share 1 allele IBD)")
    ax.set_xlim(-0.05, 1.08); ax.set_ylim(-0.05, 1.08); ax.legend(fontsize=7.5, loc="upper right", frameon=True)
    ax.set_title("--genome: every pair of animals on the Z0-Z1 plane")
    save("06_ibd_z0z1")

    # real cattle: GRM vs pedigree
    A = np.load(f"{W}/real/A_ped.npy"); G = np.loadtxt(f"{W}/real/c_grm.rel")
    iu = np.triu_indices(len(A), 1); a, gg = A[iu], G[iu]
    rng = np.random.default_rng(1); keep = (a > 0.05) | (rng.random(len(a)) < 0.03)
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.scatter(a[keep], gg[keep], s=6, alpha=0.35, color=BLUE)
    ax.plot([0, 0.5], [0, 0.5], color=RED, ls="--", label="genomic = pedigree")
    ax.set_xlabel("Pedigree relationship (from the family tree)"); ax.set_ylabel("Genomic relationship (--make-rel)")
    ax.set_title("Real cattle: half-sibs share 0.25 on average, but each pair differs")
    ax.legend(frameon=False, fontsize=9)
    save("06_real_grm_vs_pedigree")

# ---------------------------------------------------------------- 7. LD
def ld_decay(path, maxkb, step):
    d = pd.read_csv(path, sep=r"\s+", usecols=["CHR_A", "BP_A", "CHR_B", "BP_B", "R2"])
    d = d[d.CHR_A == d.CHR_B]; d["kb"] = (d.BP_B - d.BP_A) / 1000
    d = d[d.kb <= maxkb]; d["bin"] = (d.kb // step) * step + step / 2
    return d.groupby("bin").R2.mean()

def fig_ld():
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for path, lab, col, mk, st in [(f"{W}/figdata/ld_sim.ld", "simulated cattle (3 breeds)", GREY, 1000, 100),
                                   (f"{W}/real/c_ld.ld", "real cattle (one population)", BLUE, 2000, 100),
                                   (f"{W}/real/s_ld.ld", "real sheep (two breeds)", GREEN, 2000, 100)]:
        s = ld_decay(path, mk, st); ax.plot(s.index, s.values, "o-", ms=3.5, color=col, label=lab)
    ax.set_xlabel("Distance between the two SNPs (kb)"); ax.set_ylabel("Average r²")
    ax.set_title("LD decay: nearby SNPs are more correlated"); ax.legend(frameon=False)
    save("07_ld_decay")

    m = np.loadtxt(f"{W}/figdata/ld14.ld"); snps = [l.strip() for l in open(f"{W}/figdata/ld14.snplist")]
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    im = ax.imshow(np.tril(m), cmap="Reds", vmin=0, vmax=1)
    i = snps.index("snp02145") if "snp02145" in snps else None
    if i is not None:
        ax.axhline(i, color=BLUE, lw=0.8); ax.axvline(i, color=BLUE, lw=0.8)
        ax.text(i + 1.5, i - 2, "snp02145", color=BLUE, fontsize=8.5)
    ax.set_xlabel("SNPs on chromosome 14 (23.5 to 26.5 Mb)"); ax.set_ylabel("SNPs")
    plt.colorbar(im, ax=ax, fraction=0.046, label="r²"); ax.set_title("LD heatmap (--r2 square)")
    save("07_ld_heatmap")

    # schematic: mutation and recombination
    fig, ax = plt.subplots(figsize=(10, 3.4)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 34)
    def chrom(y, colours, label, mut=None):
        x = 3
        for w, c in colours:
            ax.add_patch(Rectangle((x, y), w, 3, color=c)); x += w
        ax.text(1, y + 1.5, label, ha="right", va="center", fontsize=8.5)
        if mut is not None: ax.plot(mut, y + 1.5, marker="*", ms=12, color=RED)
    ax.text(3, 31.5, "1. A new mutation (star) appears on one chromosome", fontsize=10, fontweight="bold")
    chrom(25, [(40, BLUE)], "", mut=23)
    ax.text(3, 19.5, "2. Many generations later: recombination has cut the chromosome into pieces", fontsize=10, fontweight="bold")
    chrom(13, [(12, SKY), (16, BLUE), (12, GREEN)], "", mut=23)
    chrom(7, [(8, ORANGE), (22, BLUE), (10, PURPLE)], "", mut=23)
    chrom(1, [(15, YELLOW), (13, BLUE), (12, SKY)], "", mut=23)
    ax.text(48, 9, "SNPs close to the star are still inherited with it\n(high LD); SNPs far away are not (low LD).", fontsize=10, va="center")
    save("07_ld_schematic")

# ---------------------------------------------------------------- 8. PCA
def fig_pca():
    # toy 2-SNP PCA
    rng = np.random.default_rng(3)
    a = rng.multivariate_normal([0, 0], [[1, 0.85], [0.85, 1]], 150)
    vals, vecs = np.linalg.eigh(np.cov(a.T)); order = vals.argsort()[::-1]; vals, vecs = vals[order], vecs[:, order]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    axs[0].scatter(a[:, 0], a[:, 1], s=12, color=GREY, alpha=0.7)
    for k, col, name in [(0, RED, "PC1"), (1, BLUE, "PC2")]:
        v = vecs[:, k] * np.sqrt(vals[k]) * 2
        axs[0].annotate("", xy=(v[0], v[1]), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", lw=2.5, color=col))
        axs[0].text(v[0] * 1.12, v[1] * 1.12, f"{name}\neigenvalue {vals[k]:.2f}", color=col, fontsize=9, ha="center")
    axs[0].set_aspect("equal"); axs[0].set_xlabel("SNP A (standardised genotype)"); axs[0].set_ylabel("SNP B (standardised genotype)")
    axs[0].set_title("Eigenvectors = new directions;\neigenvalues = spread along them")
    sc = a @ vecs
    axs[1].scatter(sc[:, 0], sc[:, 1], s=12, color=GREY, alpha=0.7); axs[1].set_aspect("equal")
    axs[1].set_xlabel(f"PC1 score ({vals[0]/vals.sum():.0%} of variance)"); axs[1].set_ylabel(f"PC2 score ({vals[1]/vals.sum():.0%})")
    axs[1].set_title("Same animals, rotated onto the PCs")
    save("08_pca_toy")

    ev = np.loadtxt(f"{W}/figdata/pca_raw.eigenval")
    trace = np.trace(np.loadtxt(f"{W}/figdata/rel_pr.rel"))
    fig, ax = plt.subplots(figsize=(7, 3.8))
    pct = ev / trace * 100
    ax.bar(np.arange(1, 11), pct, color=[BLUE, BLUE] + [GREY] * 8)
    for i, p in enumerate(pct): ax.text(i + 1, p + 0.15, f"{ev[i]:.1f}", ha="center", fontsize=8.5)
    ax.set_xticks(range(1, 11)); ax.set_xlabel("Principal component"); ax.set_ylabel("% of total genetic variance")
    ax.set_title("Scree plot (numbers on bars = eigenvalues; total = 384.5)")
    save("08_pca_scree")

    e = pd.read_csv(f"{W}/figdata/pca_raw.eigenvec", sep="\t"); e["breed"] = e.IID.map(breed)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, (x, y) in zip(axs, [("PC1", "PC2"), ("PC1", "PC3")]):
        for b in ["ANG", "HOL", "NDA"]:
            d = e[e.breed == b]; ax.scatter(d[x], d[y], s=14, color=BREED_COL[b], alpha=0.65, label=BREED_NAME[b])
        for iid, t in [("ANG_066", "ANG_066 (labelled Angus)"), ("HOL_077", "HOL_077 (crossbred)")]:
            r = e[e.IID == iid].iloc[0]; ax.scatter(r[x], r[y], s=80, facecolor="none", edgecolor="black", lw=1.5)
            ax.annotate(t, (r[x], r[y]), xytext=(-120 if r[x] > 0 else 8, 12), textcoords="offset points", fontsize=8.5)
        ax.axhline(0, color=GREY, ls=":"); ax.axvline(0, color=GREY, ls=":")
        ax.set_xlabel(f"{x} ({pct[int(x[2:])-1]:.1f}%)"); ax.set_ylabel(f"{y} ({pct[int(y[2:])-1]:.1f}%)")
    axs[0].legend(frameon=False, fontsize=9); axs[0].set_title("PC1 vs PC2: three breeds, three clouds"); axs[1].set_title("PC1 vs PC3: PC3 picks out a few Holstein families")
    save("08_pca_breeds")

    v = pd.read_csv(f"{W}/figdata/pca_raw.eigenvec.var", sep="\t")
    v["idx"] = np.arange(len(v)); fig, ax = plt.subplots(figsize=(10, 3.2))
    cols = [BLUE if int(c) % 2 else SKY for c in v.CHR]
    ax.scatter(v.idx, v.PC1.abs(), s=4, c=cols)
    ax.set_xlabel("SNPs in genome order (alternating colours = chromosomes)"); ax.set_ylabel("|PC1 loading|")
    ax.set_title("SNP loadings (var-wts): spread evenly along the genome, no single region drives PC1")
    save("08_pca_loadings")

    s = pd.read_csv(f"{W}/real/s_pca.eigenvec", sep="\t"); sev = np.loadtxt(f"{W}/real/s_pca.eigenval"); sp = sev / np.trace(np.loadtxt(f"{W}/real/s_rel.rel")) * 100
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for b, c in [("Jacobs", PURPLE), ("Navajo-Churro", ORANGE)]:
        d = s[s.FID == b]; ax.scatter(d.PC1, d.PC2, s=22, color=c, alpha=0.8, label=f"{b} (n={len(d)})")
    ax.axhline(0, color=GREY, ls=":"); ax.axvline(0, color=GREY, ls=":")
    ax.set_xlabel(f"PC1 ({sp[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({sp[1]:.1f}%)"); ax.legend(frameon=False)
    ax.set_title("Real sheep (Kijas et al. 2016): two breeds split on PC1")
    save("08_pca_real_sheep")

# ---------------------------------------------------------------- 9. association
def manhattan(ax, df, pcol, title, highlight=()):
    df = df.dropna(subset=[pcol]).copy(); df["lp"] = -np.log10(df[pcol])
    x = np.arange(len(df)); cols = [BLUE if int(c) % 2 else SKY for c in df.CHR]
    ax.scatter(x, df.lp, s=5, c=cols)
    for s in highlight:
        m = df.SNP == s
        if m.any(): ax.scatter(x[m.values], df.lp[m], s=50, color=RED, zorder=3); ax.annotate(s, (x[m.values][0], df.lp[m].values[0]), xytext=(5, 5), textcoords="offset points", fontsize=8.5, color=RED)
    ax.axhline(-np.log10(0.05 / len(df)), color=RED, ls="--", lw=1)
    ticks = df.groupby("CHR").apply(lambda g: x[df.index.get_indexer(g.index)].mean())
    ax.set_xticks(ticks.values[::2], [str(int(c)) for c in ticks.index[::2]], fontsize=8)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("-log10(p)"); ax.set_title(title)

def qq(ax, p, label, col):
    p = np.sort(p.dropna().values); e = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
    ax.scatter(e, -np.log10(p), s=6, color=col, label=label)

def fig_assoc():
    a = ws(f"{W}/ch3/a_basic.assoc").reset_index(drop=True); l = ws(f"{W}/ch3/lg_pc.assoc.logistic").reset_index(drop=True)
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.6))
    manhattan(axs[0], a, "P", "Naive --assoc: lambda = 10.2, hundreds of false hits", ("snp00706",))
    manhattan(axs[1], l, "P", "--logistic with PC1 and PC2: lambda = 1.04", ("snp00706",))
    plt.tight_layout(); save("09_manhattan_brd")
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    qq(ax, a.P, "naive --assoc (λ = 10.2)", RED); qq(ax, l.P, "--logistic + PCs (λ = 1.04)", BLUE)
    lim = ax.get_xlim()[1]; ax.plot([0, lim], [0, lim], color="black", lw=1)
    ax.set_xlabel("Expected -log10(p)"); ax.set_ylabel("Observed -log10(p)"); ax.legend(frameon=False, fontsize=9)
    ax.set_title("QQ plot: points should follow the line"); save("09_qq_brd")

    q = ws(f"{W}/ch4/q_basic.qassoc").reset_index(drop=True); li = ws(f"{W}/ch4/lin.assoc.linear").reset_index(drop=True)
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.6))
    manhattan(axs[0], q, "P", "Weaning weight, naive --assoc (qassoc): lambda = 20.2", ("snp02145",))
    manhattan(axs[1], li, "P", "Weaning weight, --linear with PCs, sex, feedlot, age: lambda = 1.10", ("snp02145", "snp00891"))
    plt.tight_layout(); save("09_manhattan_wwt")

def fig_traits():
    raw = ws(f"{W}/figdata/g02145.raw"); ph = ws(f"{W}/data/pheno.txt")
    d = raw.merge(ph, on=["FID", "IID"]); col = [c for c in raw.columns if c.startswith("snp02145")][0]
    d = d[(d.WWT > 0)]; d["breed"] = d.IID.map(breed)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for k, b in enumerate(["ANG", "HOL", "NDA"]):
        means = [d[(d.breed == b) & (d[col] == gt)].WWT.mean() for gt in (0, 1, 2)]
        ax.plot([0, 1, 2], means, "o-", color=BREED_COL[b], label=BREED_NAME[b], lw=2)
    ax.set_xticks([0, 1, 2], [f"0 copies of {col.split('_')[1]}", "1 copy", "2 copies"])
    ax.set_ylabel("Mean weaning weight (kg)"); ax.legend(frameon=False)
    ax.set_title("snp02145: same slope in every breed, very different starting levels")
    save("10_wwt_by_genotype")

    g = ws(f"{W}/figdata/g01540.raw"); cv = ws(f"{W}/data/covar.txt")
    d = g.merge(ph, on=["FID", "IID"]).merge(cv, on=["FID", "IID"]); gc = [c for c in g.columns if c.startswith("snp01540")][0]
    d = d[d.WWT > 0]; d["breed"] = d.IID.map(breed)
    d["resid"] = d.WWT - d.groupby("breed").WWT.transform("mean")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for env, c, lab in [(0, GREEN, "pasture"), (1, ORANGE, "feedlot")]:
        m = [d[(d.FEEDLOT == env) & (d[gc] == k)].resid.mean() for k in (0, 1, 2)]
        ax.plot([0, 1, 2], m, "o-", color=c, lw=2, label=lab)
    ax.set_xticks([0, 1, 2], [f"0 x {gc.split('_')[1]}", "1 copy", "2 copies"]); ax.set_ylabel("Weaning weight, deviation from breed mean (kg)")
    ax.legend(frameon=False); ax.set_title("Gene x environment: snp01540 matters only in the feedlot")
    save("10_gxe_lines")

    # forest plot meta-analysis
    rows = []
    for b in ["ANG", "HOL", "NDA"]:
        x = ws(f"{W}/ch7/mk_{b}.assoc.linear"); r = x[x.SNP == "snp02145"].iloc[0]; rows.append((BREED_NAME[b], r.BETA, r.L95, r.U95))
    ses = [(r[1], (r[3] - r[2]) / (2 * 1.96)) for r in rows]
    w = [1 / s ** 2 for _, s in ses]; b = sum(wi * bi for wi, (bi, _) in zip(w, ses)) / sum(w); se = (1 / sum(w)) ** 0.5
    rows.append(("Meta-analysis (fixed)", b, b - 1.96 * se, b + 1.96 * se))
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    for i, (n, b, lo, hi) in enumerate(rows[::-1]):
        diamond = n.startswith("Meta")
        ax.plot([lo, hi], [i, i], color="black"); ax.plot(b, i, "D" if diamond else "s", ms=10 if diamond else 7, color=RED if diamond else BLUE)
        ax.text(hi + 0.5, i, f"{b:.1f} kg", va="center", fontsize=9)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows[::-1]]); ax.axvline(0, color=GREY, ls=":")
    ax.set_xlabel("Effect of each T allele on weaning weight (kg, 95% CI)"); ax.set_title("Forest plot: snp02145 in three breeds")
    save("10_forest_meta")

def fig_epi_cnv_roh():
    # epistasis: proportion of PKE cases by genotype pair
    ph = ws(f"{W}/ch7/pke.txt"); g = ws(f"{W}/ch7/epi_raw.raw")
    d = g.merge(ph, on=["FID", "IID"]); c1, c2 = [c for c in g.columns if c.startswith("snp")]
    tab = d.groupby([c1, c2]).PKE.apply(lambda s: (s == 2).mean()).unstack()
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    im = ax.imshow(tab.values, cmap="RdBu_r", vmin=0, vmax=1)
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            n = ((d[c1] == tab.index[i]) & (d[c2] == tab.columns[j])).sum()
            ax.text(j, i, f"{tab.values[i, j]:.0%}\n(n={n})", ha="center", va="center", fontsize=9)
    ax.set_xticks(range(3), [f"{int(v)}" for v in tab.columns]); ax.set_yticks(range(3), [f"{int(v)}" for v in tab.index])
    ax.set_xlabel(f"copies of {c2.split('_')[1]} at snp03388"); ax.set_ylabel(f"copies of {c1.split('_')[1]} at snp01294")
    plt.colorbar(im, fraction=0.046, label="share of cases"); ax.set_title("Epistasis: risk is high only in two\nopposite corners")
    save("11_epistasis_grid")

    s = ws(f"{W}/ch6/loci.cnv.summary"); s = s[s.CHR == 7]
    fig, ax = plt.subplots(figsize=(8, 3.4))
    ax.step(s.BP / 1e6, s.AFF, where="post", color=RED, label="cases"); ax.step(s.BP / 1e6, s.UNAFF, where="post", color=BLUE, label="controls")
    ax.set_xlabel("Chromosome 7 position (Mb)"); ax.set_ylabel("Animals with a CNV here"); ax.legend(frameon=False)
    ax.set_title("CNV counts along chromosome 7: one deletion hotspot near 40 Mb")
    save("12_cnv_chr7")

    L = 2507000
    sim = ws(f"{W}/ch2/rohC.hom.indiv") if os.path.exists(f"{W}/ch2/rohC.hom.indiv") else None
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    if sim is not None:
        sim["FROH"] = sim.KB / L; sim["breed"] = sim.IID.map(breed)
        for b in ["ANG", "HOL", "NDA"]:
            axs[0].hist(sim[sim.breed == b].FROH, bins=np.arange(0, 0.32, 0.01), color=BREED_COL[b], alpha=0.7, label=BREED_NAME[b])
        v = sim[sim.IID == "NDA_040"].FROH.values[0]
        axs[0].annotate("NDA_040", (v, 1), xytext=(v - 0.08, 40), arrowprops=dict(arrowstyle="->"), fontsize=9)
        axs[0].legend(frameon=False, fontsize=9); axs[0].set_xlabel("F_ROH"); axs[0].set_ylabel("Animals")
        axs[0].set_title("Simulated cattle (ROH ≥ 10 Mb, ≥ 30 SNPs)")
    r = ws(f"{W}/real/s_roh.hom.indiv"); bim = ws(f"{W}/real/sheep_qc.bim", header=None)
    Ls = bim.groupby(0)[3].max().sum() / 1000
    r["FROH"] = r.KB / Ls
    for b, c in [("Jacobs", PURPLE), ("Navajo-Churro", ORANGE)]:
        axs[1].hist(r[r.FID == b].FROH, bins=np.arange(0, 0.6, 0.02), color=c, alpha=0.75, label=b)
    axs[1].legend(frameon=False, fontsize=9); axs[1].set_xlabel("F_ROH (chromosomes 2 and 24 only)"); axs[1].set_ylabel("Animals")
    axs[1].set_title("Real sheep (ROH ≥ 1 Mb, ≥ 30 SNPs)")
    save("13_froh")

def fig_pca_proj():
    e = pd.read_csv(f"{W}/figdata/pca_raw.eigenvec", sep="\t")
    p = pd.read_csv(f"{W}/figdata/proj.sscore", sep="\t").rename(columns={"#FID": "FID"})
    fam = ws(f"{W}/ch1/cattle.fam", header=None, names=["FID", "IID", "PAT", "MAT", "SEX", "PHE"])
    calves = set(fam[fam.PAT != "0"].IID); parents = set(fam[fam.FID.str.contains("_F")].IID)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, df, x, y, title in [(axs[0], e, "PC1", "PC3", "PCs computed on all 361 animals"),
                                (axs[1], p, "PC1_AVG", "PC3_AVG", "PCs computed on 301 founders, everyone projected")]:
        fams = df.IID.isin(parents)
        ax.scatter(df[~fams][x], df[~fams][y], s=12, color=GREY, alpha=0.5, label="unrelated animals")
        ax.scatter(df[fams][x], df[fams][y], s=18, color=ORANGE, alpha=0.9, label="members of the 30 families")
        ax.set_xlabel("PC1"); ax.set_ylabel("PC3"); ax.set_title(title)
    axs[0].legend(frameon=False, fontsize=9)
    save("08_pca_projection")

def fig_params():
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    g = [("0.01", 967, 25), ("0.02", 49, 25), ("0.05", 33, 25), ("0.10", 18, 18), ("0.20", 5, 5)]
    x = range(len(g))
    axs[0].bar(x, [a for _, a, _ in g], color=GREY, label="good SNPs removed too")
    axs[0].bar(x, [b for _, _, b in g], color=ORANGE, label="planted bad SNPs")
    axs[0].set_xticks(list(x), [k for k, _, _ in g]); axs[0].set_xlabel("--geno threshold")
    axs[0].set_ylabel("SNPs removed (of 4,800)"); axs[0].set_title("--geno: stricter is not always better")
    for i, (_, a, _) in enumerate(g): axs[0].text(i, a + 15, str(a), ha="center", fontsize=9)
    axs[0].legend(frameon=False, fontsize=9)
    m = [("0.01", 314), ("0.05", 802), ("0.10", 1329)]
    axs[1].bar(range(3), [v for _, v in m], color=BLUE)
    axs[1].set_xticks(range(3), [k for k, _ in m]); axs[1].set_xlabel("--maf threshold"); axs[1].set_title("--maf: SNPs removed")
    for i, (_, v) in enumerate(m): axs[1].text(i, v + 20, f"{v:,}", ha="center", fontsize=9)
    h = [("pooled\ncontrols", 63), ("pooled\nall", 114), ("pooled\n+calves", 81), ("within\nbreed", 12)]
    axs[2].bar(range(4), [v for _, v in h], color=[GREY, GREY, GREY, GREEN])
    axs[2].axhline(12, ls=":", color=ORANGE); axs[2].text(1.5, 16, "12 truly bad SNPs", color=ORANGE, fontsize=9, ha="center")
    axs[2].set_xticks(range(4), [k for k, _ in h], fontsize=10); axs[2].set_title("--hwe 1e-6: SNPs removed")
    for i, (_, v) in enumerate(h): axs[2].text(i, v + 2, str(v), ha="center", fontsize=9)
    save("14_param_effects")

def fig_real():
    het = ws(f"{W}/real/c_het.het"); m, sd = het.F.mean(), het.F.std()
    e = ws(f"{W}/real/c_pca.eigenvec"); ev = np.loadtxt(f"{W}/real/c_pca.eigenval"); tr = np.trace(np.loadtxt(f"{W}/real/c_grm.rel"))
    hi = ws(f"{W}/real/c_roh.hom.indiv"); bim = ws(f"{W}/real/cattle_real_qc.bim", header=None)
    L = (bim.groupby(0)[3].max() - bim.groupby(0)[3].min()).sum() / 1000
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    axs[0].hist(het.F, bins=40, color=BLUE)
    for k in (-3, 3): axs[0].axvline(m + k * sd, color=RED, ls="--")
    axs[0].set_xlabel("Inbreeding coefficient F (--het)"); axs[0].set_ylabel("Animals")
    axs[0].set_title(f"F: {(abs(het.F - m) > 3 * sd).sum()} animals beyond 3 SD")
    axs[1].scatter(e.PC1, e.PC2, s=10, color=BLUE, alpha=0.6)
    axs[1].set_xlabel(f"PC1 ({ev[0] / tr * 100:.1f}%)"); axs[1].set_ylabel(f"PC2 ({ev[1] / tr * 100:.1f}%)")
    axs[1].set_title("PCA: one population, no breed clusters")
    axs[2].hist(hi.KB / L, bins=30, color=GREEN)
    axs[2].set_xlabel("F_ROH (runs of at least 2 Mb)"); axs[2].set_ylabel("Animals"); axs[2].set_title(f"F_ROH: mean {(hi.KB / L).mean():.3f}")
    save("15_real_cattle")

if __name__ == "__main__":
    for f in [fig_formats, fig_toy_missing, fig_missing, fig_maf_hwe, fig_sex_het, fig_ibd, fig_ld, fig_pca, fig_assoc, fig_traits, fig_epi_cnv_roh, fig_pca_proj, fig_params, fig_real]:
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            import traceback; print("FAIL", f.__name__, e); traceback.print_exc()
