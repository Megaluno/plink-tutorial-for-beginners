"""Simulate a second binary trait (PKE) with a pure SNP x SNP interaction
(snp01294 on chr 8 x snp03388 on chr 22; both have similar MAF in all breeds)."""
import random
random.seed(42)
rows = [l.split() for l in open('epi_raw.raw')][1:]
def num(v):
    return None if v == 'NA' else float(v)
vals = [(r[0], r[1], num(r[6]), num(r[7])) for r in rows]
ma = sum(v[2] for v in vals if v[2] is not None) / sum(1 for v in vals if v[2] is not None)
mb = sum(v[3] for v in vals if v[3] is not None) / sum(1 for v in vals if v[3] is not None)
lia = []
for f, i, a, b in vals:
    a = ma if a is None else a
    b = mb if b is None else b
    lia.append((f, i, 1.6 * (a - ma) * (b - mb) + random.gauss(0, 1)))
th = sorted(l[2] for l in lia)[int(0.6 * len(lia))]
with open('pke.txt', 'w') as o:
    o.write('FID IID PKE\n')
    for f, i, l in lia:
        o.write(f"{f} {i} {2 if l > th else 1}\n")
print('cases', sum(1 for l in lia if l[2] > th), 'of', len(lia))
