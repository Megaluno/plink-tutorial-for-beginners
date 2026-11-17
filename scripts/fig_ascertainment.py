"""Figure for the ascertainment-bias section, from the files of scripts/ascertainment_pipeline.sh."""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = sys.argv[1] if len(sys.argv) > 1 else "work_ascertainment"
OUT = sys.argv[2] if len(sys.argv) > 2 else "images"
BLUE, ORANGE, GREEN, RED, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#7f7f7f"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlesize": 12, "savefig.dpi": 150, "savefig.bbox": "tight"})
BREEDS = ["EUR1", "EUR2", "AFR"]
SETS = [("seq", "Sequence (random SNPs)", GREY), ("chipEUR", "Chip discovered in EUR1", RED),
        ("chipALL_n", "Chip discovered in all breeds", BLUE)]

def strat(s):
    return pd.read_csv(f"{W}/f_{s}.frq.strat", sep=r"\s+")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [1, 1.25]})
x = np.arange(3); w = 0.26
for k, (s, lab, col) in enumerate(SETS):
    d = strat(s)
    he = [np.mean(2 * d[d.CLST == b].MAF * (1 - d[d.CLST == b].MAF)) for b in BREEDS]
    a1.bar(x + (k - 1) * w, he, w, color=col, label=lab)
    for xi, v in zip(x + (k - 1) * w, he):
        a1.text(xi, v + 0.004, f"{v:.2f}", ha="center", fontsize=8.5)
a1.set_xticks(x); a1.set_xticklabels(BREEDS); a1.set_ylabel("Expected heterozygosity (mean 2pq)")
a1.set_title("Diversity of each breed, by SNP set"); a1.legend(frameon=False, fontsize=9, loc="upper right")
a1.set_ylim(0, 0.39)

bins = np.linspace(0, 0.5, 11)
for s, lab, col, ls in (("seq", "Sequence", GREY, "-"), ("chipEUR", "Chip discovered in EUR1", RED, "-")):
    for b, mk in (("EUR1", "o"), ("AFR", "s")):
        d = strat(s); p = d[d.CLST == b].MAF.values; m = np.minimum(p, 1 - p)
        poly = m[m > 0]
        h, _ = np.histogram(poly, bins=bins)
        a2.plot((bins[:-1] + bins[1:]) / 2, h / len(m), marker=mk, color=col, ls="-" if b == "EUR1" else "--",
                 label=f"{lab}, {b} ({np.mean(m == 0) * 100:.0f}% monomorphic)")
a2.set_xlabel("Minor allele frequency within the breed (polymorphic SNPs)"); a2.set_ylabel("Share of all SNPs in the set")
a2.set_title("Allele frequency spectrum"); a2.legend(frameon=False, fontsize=8.5)
plt.tight_layout(); plt.savefig(f"{OUT}/18_ascertainment.png"); plt.close()
print("written", f"{OUT}/18_ascertainment.png")
