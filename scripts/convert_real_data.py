"""Convert the real cattle data used in Chapter 18 to PLINK .ped/.map format.

Source: synbreedData R package (GPL-2), file data/cattle.RData
        https://github.com/cran/synbreedData
The real sheep data (Kijas et al. 2016) already come as PLINK .ped/.map files in the
detectRUNS package (GPL-3): https://github.com/cran/detectRUNS, folder inst/extdata.

Usage:  pip install rdata numpy pandas
        git clone --depth 1 https://github.com/cran/synbreedData.git
        python convert_real_data.py synbreedData/data/cattle.RData out_folder
Then:   plink --file out_folder/cattle_real --cow --make-bed --out cattle_real
"""
import sys, os, warnings
import numpy as np, rdata

warnings.filterwarnings("ignore")
src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
cat = rdata.conversion.convert(rdata.parser.parse_file(src))["cattle"]

mp = cat["map"]                                   # rownames = SNP names; chr, pos (Mb)
snps = list(mp.index)
ids = list(cat["covar"].query("genotyped").id)    # genotyped animals, in matrix order
geno = np.array(cat["geno"], dtype=object).reshape(len(ids), len(snps), order="F")  # R stores by column

with open(f"{out}/cattle_real.ped", "w") as f:
    for i, iid in enumerate(ids):
        calls = ["0 0" if g is None else f"{g[0]} {g[1]}" for g in geno[i]]
        f.write(f"CATTLE {iid} 0 0 1 -9 " + " ".join(calls) + "\n")
with open(f"{out}/cattle_real.map", "w") as f:
    for s in snps:
        bp = int(round(mp.loc[s, "pos"] * 1e6)) + 1   # Mb -> bp, starting at 1
        f.write(f"{int(mp.loc[s, 'chr'])}\t{s}\t0\t{bp}\n")

cat["pedigree"].to_csv(f"{out}/cattle_real_pedigree.csv", index=False)
print(f"{len(ids)} animals and {len(snps)} SNPs written to {out}/cattle_real.ped/.map")
