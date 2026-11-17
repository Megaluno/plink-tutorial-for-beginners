"""Figures for the phasing and imputation part of Chapter 14, from the files written by
scripts/chapters/ch14_phasing.sh (folder given as the first argument, default work/phase)."""
import sys, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = sys.argv[1] if len(sys.argv) > 1 else "work/phase"
OUT = sys.argv[2] if len(sys.argv) > 2 else "images"
BLUE, ORANGE, GREEN, RED, PURPLE, SKY, YELLOW, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#7f7f7f"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlesize": 12, "figure.dpi": 110,
                     "savefig.dpi": 150, "savefig.bbox": "tight"})

# ------------------------------------------------ 1. copying mosaic (Li-Stephens Viterbi paths)
def mosaic():
    fig, axes = plt.subplots(2, 1, figsize=(10, 3.6), sharex=True)
    pal = [BLUE, ORANGE, GREEN, RED, PURPLE, SKY, YELLOW, "#999999", "#000000", "#8c564b"]
    for ax, (fn, title) in zip(axes, [("paths_refA10.txt", "Reference of 10 animals (20 haplotypes)"),
                                      ("paths_refA.txt", "Reference of 150 animals (300 haplotypes)")]):
        d = pd.read_csv(f"{W}/{fn}", sep="\t")
        d = d[d.IID == "T_091"]
        x = d.POS.values / 1e6
        edges = np.r_[x[0] - 0.02, (x[1:] + x[:-1]) / 2, x[-1] + 0.02]
        for row, col in [(1, "HAP1_FROM"), (0, "HAP2_FROM")]:
            src = d[col].values
            names = list(dict.fromkeys(src))
            start = 0
            for k in range(1, len(src) + 1):
                if k == len(src) or src[k] != src[start]:
                    c = pal[names.index(src[start]) % len(pal)]
                    ax.barh(row, edges[k] - edges[start], left=edges[start], height=0.7, color=c, edgecolor="white", lw=0.6)
                    start = k
        ax.set_yticks([0, 1]); ax.set_yticklabels(["copy 2", "copy 1"])
        ax.set_title(title, loc="left", fontsize=11)
        ax.spines["left"].set_visible(False)
    axes[1].set_xlabel("Position on chromosome 14 (Mb)")
    fig.suptitle("Animal T_091: which reference haplotype each piece of its two chromosomes was copied from",
                 x=0.01, ha="left", fontsize=11.5, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT}/16_ls_mosaic.png"); plt.close()

# ------------------------------------------------ 2. accuracy by MAF for five reference panels
def accuracy():
    runs = [("ref_A_hd", "Same breed, 150 animals", BLUE, "-", "o"),
            ("refAB", "Same breed 150 + other breed 150", SKY, "--", "s"),
            ("refA30", "Same breed, 30 animals", GREEN, "-", "^"),
            ("refA10", "Same breed, 10 animals", ORANGE, "-", "v"),
            ("ref_B_hd", "Other breed only, 150 animals", RED, "-", "D")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for key, lab, col, ls, mk in runs:
        rows = []
        for line in open(f"{W}/eval_{key}.txt"):
            t = line.split()
            if len(t) == 5 and re.match(r"^\d\.\d\d-\d\.\d\d$", t[0]):
                rows.append((t[0], float(t[2]), float(t[3])))
        cls = [r[0] for r in rows]
        a1.plot(cls, [r[2] for r in rows], ls=ls, marker=mk, color=col, label=lab, lw=2)
        a2.plot(cls, [r[1] for r in rows], ls=ls, marker=mk, color=col, label=lab, lw=2)
    a1.set_ylabel("Mean r² (dosage vs true genotype)"); a1.set_title("Imputation accuracy, r²")
    a2.set_ylabel("Share of genotypes called correctly"); a2.set_title("Concordance (can mislead)")
    for a in (a1, a2):
        a.set_xlabel("Minor allele frequency of the imputed SNP"); a.set_ylim(0, 1.03)
    a1.legend(fontsize=8.5, loc="lower right", frameon=False)
    plt.tight_layout(); plt.savefig(f"{OUT}/16_imputation_accuracy.png"); plt.close()

# ------------------------------------------------ 3. association before and after imputation
def manhattan():
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.2), sharex=True, sharey=True)
    qtl = 24920712
    for ax, (fn, title, col) in zip(axes, [("demo/chip.WWT.glm.linear", "Chip SNPs only (294)", GREY),
                                           ("demo/imp.WWT.glm.linear", "Chip + imputed SNPs with DR2 ≥ 0.8 (2,050)", BLUE)]):
        d = pd.read_csv(f"{W}/{fn}", sep="\t")
        ax.scatter(d.POS / 1e6, -np.log10(d.P), s=9, color=col)
        ax.axvline(qtl / 1e6, color=RED, ls=":", lw=1.4)
        ax.axhline(-np.log10(0.05 / 2050), color="#333", ls="--", lw=1)
        top = d.loc[d.P.idxmin()]
        ax.annotate(f"{top.ID}  p = {top.P:.1e}", (top.POS / 1e6, -np.log10(top.P)), xytext=(8, -2),
                    textcoords="offset points", fontsize=9)
        ax.set_title(title, loc="left", fontsize=11); ax.set_ylabel("-log10(p)")
    axes[0].text(qtl / 1e6 + 0.08, 5.3, "true QTL", color=RED, fontsize=9)
    axes[1].set_xlabel("Position on chromosome 14 (Mb)")
    plt.tight_layout(); plt.savefig(f"{OUT}/16_imputed_gwas.png"); plt.close()

# ------------------------------------------------ 4. DR2 against true r2
def dr2():
    import gzip
    def read(fn):
        out = {}
        ids = None
        for l in gzip.open(fn, "rt"):
            if l.startswith("##"): continue
            t = l.rstrip("\n").split("\t")
            if l.startswith("#"): ids = t[9:]; continue
            out[t[2]] = t
        return ids, out
    ii, imp = read(f"{W}/imp_refA10.vcf.gz"); ti, tru = read(f"{W}/truth_target_hd.vcf.gz")
    col = [ti.index(i) for i in ii]
    x, y = [], []
    for sid, t in imp.items():
        if "IMP" not in t[7]: continue
        ds = np.array([float(c.split(":")[1]) for c in t[9:]])
        g = np.array([int(tru[sid][9 + k][0]) + int(tru[sid][9 + k][2]) for k in col])
        if g.std() == 0 or ds.std() == 0: continue
        info = dict(kv.split("=") for kv in t[7].split(";") if "=" in kv)
        x.append(float(info["DR2"])); y.append(np.corrcoef(ds, g)[0, 1] ** 2)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(x, y, s=8, alpha=0.5, color=BLUE)
    ax.plot([0, 1], [0, 1], color=GREY, lw=1)
    ax.axvline(0.8, color=RED, ls="--", lw=1)
    ax.set_xlabel("DR2 reported by the imputation"); ax.set_ylabel("True r² (needs the real genotypes)")
    ax.set_title("Reference of 10 animals")
    plt.tight_layout(); plt.savefig(f"{OUT}/16_dr2_vs_true.png"); plt.close()

if __name__ == "__main__":
    mosaic(); accuracy(); manhattan(); dr2()
    print("figures written to", OUT)
