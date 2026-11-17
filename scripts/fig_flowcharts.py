"""Workflow diagrams for the book, drawn as figures (replacing the earlier text flowcharts)."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

OUT = sys.argv[1] if len(sys.argv) > 1 else "images"
BLUE, ORANGE, GREEN, RED, PURPLE, SKY, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#7f7f7f"
plt.rcParams.update({"font.size": 10, "savefig.dpi": 150, "savefig.bbox": "tight", "figure.dpi": 110})

def canvas(w, h, xmax, ymax):
    fig, ax = plt.subplots(figsize=(w, h)); ax.axis("off"); ax.set_xlim(0, xmax); ax.set_ylim(-0.5, ymax)
    return fig, ax

def box(ax, x, y, w, h, title, lines=(), col=BLUE, mono=True, tsize=10.5, lsize=9):
    """x, y = centre of the box"""
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                fc=col, ec="none", alpha=0.16))
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                fc="none", ec=col, lw=1.6))
    n = len(lines)
    top = y + (n * 1.0) / 2 + (0.45 if n else 0)
    ax.text(x, top if n else y, title, ha="center", va="center", fontweight="bold", fontsize=tsize)
    for i, t in enumerate(lines):
        ax.text(x, top - 1.15 - i * 1.0, t, ha="center", va="center", fontsize=lsize,
                family="monospace" if mono else None)

def diamond(ax, x, y, w, h, text, col=PURPLE):
    ax.add_patch(Polygon([(x, y + h / 2), (x + w / 2, y), (x, y - h / 2), (x - w / 2, y)], closed=True,
                         fc=col, alpha=0.16, ec="none"))
    ax.add_patch(Polygon([(x, y + h / 2), (x + w / 2, y), (x, y - h / 2), (x - w / 2, y)], closed=True,
                         fc="none", ec=col, lw=1.6))
    ax.text(x, y, text, ha="center", va="center", fontsize=9.5, fontweight="bold")

def arrow(ax, x1, y1, x2, y2, label=None, lx=0, ly=0, col="#333", rad=0.0):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color=col, connectionstyle=f"arc3,rad={rad}"))
    if label:
        ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, ha="center", va="center", fontsize=9,
                style="italic", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png"); plt.close(fig)

# 1. merge workflow (Chapter 5)
def merge():
    fig, ax = canvas(11.5, 4.6, 47, 18)
    box(ax, 5, 12, 8.5, 3.6, "Two filesets", ["cattle", "batch2"], GREY)
    box(ax, 17, 12, 9, 3.6, "Merge", ["--bmerge batch2"], BLUE)
    box(ax, 17, 3.5, 11, 4.2, "Strand problem", ["Error: 3+ alleles", "list: *-merge.missnp"], RED)
    box(ax, 31, 3.5, 10, 3.6, "Flip file 2", ["--flip *.missnp"], ORANGE)
    box(ax, 31, 12, 9.5, 3.6, "Merged file", ["merged.bed/.bim/.fam"], GREEN)
    box(ax, 42.5, 12, 7.4, 5.6, "Check", ["remove A/T,", "C/G SNPs;", "--geno;", "--flip-scan"], GREEN, lsize=8.5)
    arrow(ax, 9.6, 12, 12.3, 12)
    arrow(ax, 21.7, 12, 26, 12, "works", ly=0.9)
    arrow(ax, 17, 10, 17, 5.8, "fails", lx=1.5)
    arrow(ax, 22.7, 3.5, 25.8, 3.5)
    arrow(ax, 31, 5.5, 20.5, 10.1, "merge again", lx=1.5, ly=-0.9, rad=0.25)
    arrow(ax, 35.9, 12, 38.5, 12)
    save(fig, "flow_merge")

# 2. two-round missingness filter (Chapter 6)
def missingness():
    fig, ax = canvas(12, 3.2, 56, 12)
    steps = [("Raw data", ["361 animals", "4,800 SNPs"], GREY),
             ("Round 1, SNPs", ["--geno 0.10"], ORANGE),
             ("Round 1, animals", ["--mind 0.10"], ORANGE),
             ("Round 2, animals", ["--mind 0.03"], BLUE),
             ("Round 2, SNPs", ["--geno 0.01"], BLUE),
             ("Clean data", ["355 animals", "4,503 SNPs"], GREEN)]
    xs = [4, 13.5, 23.5, 33.5, 43.5, 52.5]
    for (t, l, c), x in zip(steps, xs):
        box(ax, x, 6, 7.4 if x not in (4, 52.5) else 6.2, 4.2, t, l, c, tsize=9.5)
    half = {x: (3.1 if x in (4, 52.5) else 3.7) + 0.35 for x in xs}
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, a + half[a], 6, b - half[b], 6)
    ax.text(18.5, 10.6, "lenient: only obvious failures", ha="center", fontsize=9, style="italic", color="#8a5a00")
    ax.text(38.5, 10.6, "strict: animals before SNPs", ha="center", fontsize=9, style="italic", color="#004f7c")
    ax.text(13.5, 1.6, "-18 SNPs", ha="center", fontsize=9); ax.text(23.5, 1.6, "-4 animals", ha="center", fontsize=9)
    ax.text(33.5, 1.6, "-2 animals", ha="center", fontsize=9); ax.text(43.5, 1.6, "-279 SNPs", ha="center", fontsize=9)
    save(fig, "flow_missingness")

# 3. choosing a case/control test (Chapter 12)
def assoc():
    fig, ax = canvas(11, 6.4, 46, 26)
    box(ax, 23, 23.5, 10, 2.6, "Case/control trait", [], GREY)
    diamond(ax, 23, 18.3, 14, 5, "More than one breed,\nherd or batch?")
    diamond(ax, 9, 11.5, 12, 4.6, "Rare SNPs or\nsmall study?")
    box(ax, 4.5, 3, 8.4, 3.4, "Exact test", ["--assoc fisher"], ORANGE)
    box(ax, 14.5, 3, 10, 4.2, "Allelic / genotype tests", ["--assoc, --model", "+ --adjust (check λ)"], ORANGE, tsize=9.5)
    box(ax, 35, 11.5, 15, 4.2, "Model the structure", ["--logistic --covar pcs.txt", "or --mh --within breed.clst"], BLUE)
    box(ax, 29, 3, 11, 3.6, "Same effect everywhere?", ["--bd, --homog"], GREEN, tsize=9.5)
    box(ax, 41, 3, 9.5, 3.6, "Additive?", ["--logistic --genotypic"], GREEN, tsize=9.5)
    arrow(ax, 23, 22.1, 23, 20.9)
    arrow(ax, 16.2, 17.3, 11.5, 13.5, "no", lx=-1.3)
    arrow(ax, 29.8, 17.3, 33, 13.7, "yes", lx=1.3)
    arrow(ax, 6.4, 9.6, 4.8, 4.8, "yes", lx=-1.3)
    arrow(ax, 11.6, 9.6, 13.6, 5.2, "no", lx=1.3)
    arrow(ax, 32, 9.3, 29.5, 4.9)
    arrow(ax, 38, 9.3, 40.5, 4.9)
    save(fig, "flow_assoc_choice")

# 4. imputation idea (Chapter 14)
def imputation():
    fig, ax = canvas(11, 4.4, 46, 17)
    def hap(x, y, alleles, hi=None, col="#333"):
        for i, a in enumerate(alleles):
            c = RED if i == hi else col
            ax.text(x + i * 2.0, y, a, ha="center", va="center", family="monospace", fontsize=13,
                    fontweight="bold" if i == hi else "normal", color=c)
            if i < len(alleles) - 1: ax.text(x + i * 2.0 + 1.0, y, "-", ha="center", va="center", fontsize=12, color=GREY)
    ax.text(2, 15.5, "Reference panel (dense, phased)", fontweight="bold", fontsize=10.5, color="#004f7c")
    ax.text(28, 15.5, "Study animal (chip only)", fontweight="bold", fontsize=10.5, color="#8a5a00")
    for k, (al, y) in enumerate([("AGTC", 12.5), ("CTGA", 10.0), ("AGTA", 7.5)]):
        ax.text(1.5, y, f"haplotype {k + 1}", fontsize=9, va="center")
        hap(10, y, list(al), hi=2)
    ax.add_patch(FancyBboxPatch((9, 11.5), 7.9, 2, boxstyle="round,pad=0.2", fc=GREEN, alpha=0.18, ec=GREEN))
    ax.text(32, 12.5, "copy 1", fontsize=9, va="center", ha="right"); hap(34, 12.5, list("AG?C"), hi=2, col="#333")
    ax.text(32, 9.0, "copy 2", fontsize=9, va="center", ha="right"); hap(34, 9.0, list("CT?A"), hi=2, col="#333")
    arrow(ax, 27.5, 12.5, 17.8, 12.5, "copy 1 matches haplotype 1", ly=1.0)
    arrow(ax, 27.5, 9.0, 17.8, 9.9, "copy 2 matches haplotype 2", ly=-1.1)
    box(ax, 23, 3, 30, 3.6, "Imputed genotype at SNP 3: T / G", ["dosage of T ≈ 1.0, with a probability for each genotype"], GREEN, mono=False)
    save(fig, "flow_imputation")

# 5. order of filters inside one PLINK 1.9 command (Chapter 17)
def order():
    fig, ax = canvas(11.5, 2.8, 52, 10)
    steps = [("1", ["--keep / --remove", "--extract / --exclude"]), ("2", ["--mind"]), ("3", ["--geno"]),
             ("4", ["--hwe"]), ("5", ["--maf"]), ("6", ["analysis or", "--make-bed"])]
    xs = [5.5, 15.5, 23.5, 31.5, 39.5, 47.5]
    ws = [10, 5.6, 5.6, 5.6, 5.6, 7.6]
    for (n, l), x, w in zip(steps, xs, ws):
        box(ax, x, 5, w, 4, f"step {n}", l, BLUE if n != "6" else GREEN, tsize=9.5, lsize=8.8)
    for i in range(5):
        arrow(ax, xs[i] + ws[i] / 2 + 0.3, 5, xs[i + 1] - ws[i + 1] / 2 - 0.3, 5)
    ax.text(26, 9.3, "whatever order you type the flags in", ha="center", fontsize=9.5, style="italic")
    save(fig, "flow_filter_order")

# 6. the QC pipeline (Chapter 19)
def pipeline():
    fig, ax = canvas(10.5, 11, 40, 46.5)
    rows = [("Raw data", ["361 animals, 4,800 SNPs"], GREY, 43.5, 3.0),
            ("0. Look (nothing removed yet)", ["--missing  --freq  --hardy (within breed)", "--test-missing  --het  --check-sex",
                                               "--mendel  --genome  --pca"], ORANGE, 36.3, 5.0),
            ("1. Decide: write the lists", ["remove_animals.txt  sex_update.txt", "snp_exclude.txt  breed_fixed.clst"], PURPLE, 28.6, 4.0),
            ("2a. Animals first", ["--remove remove_animals.txt --update-sex sex_update.txt"], BLUE, 22.3, 3.0),
            ("2b. Then SNPs", ["--exclude snp_exclude.txt --geno 0.05"], BLUE, 17.1, 3.0),
            ("2c. Safety check", ["--mind 0.05   (expect 0 removed)"], BLUE, 11.9, 3.0),
            ("2d. Rare SNPs", ["--maf 0.01"], BLUE, 6.7, 3.0),
            ("Clean data: cattle_qc", ["351 animals, 4,443 SNPs  ->  prune, PCA for association"], GREEN, 1.5, 3.0)]
    for t, l, c, y, h in rows:
        box(ax, 20, y, 30, h, t, l, c, tsize=10, lsize=8.8)
    for (_, _, _, y1, h1), (_, _, _, y2, h2) in zip(rows[:-1], rows[1:]):
        arrow(ax, 20, y1 - h1 / 2 - 0.35, 20, y2 + h2 / 2 + 0.35)
    save(fig, "flow_pipeline")

if __name__ == "__main__":
    for f in (merge, missingness, assoc, imputation, order, pipeline):
        f(); print("ok", f.__name__)
