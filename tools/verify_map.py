"""Verification: render the new map's views headless for visual comparison
against the old map, and run determinism / statistics checks.

Run with:  /Users/ashwinnimmal/opt/anaconda3/bin/python tools/verify_map.py
"""
import collections
import os
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "data")

import starmap2
from starmap2 import fig, galaxy, show_galaxy_map, show_sector_map, show_system_view
from stellarforge import CLASS_NAMES, Galaxy


def save(name):
    fig.savefig(os.path.join(OUT, name), dpi=100, facecolor="#050510")
    print("wrote data/" + name)


# --- visual renders -----------------------------------------------------
show_galaxy_map()
save("verify_galaxy.png")

sid = int(np.argmax(galaxy.density))          # a busy arm sector
show_sector_map(sid)
save("verify_sector.png")

stars = galaxy.generate_sector_stars(sid)
star_idx = int(np.flatnonzero(stars["habitable"])[0])
show_system_view(sid, star_idx, stars)
save("verify_system.png")

# --- determinism --------------------------------------------------------
g2 = Galaxy()
s1 = galaxy.generate_sector_stars(sid)
s2 = g2.generate_sector_stars(sid)
assert all(np.array_equal(s1[k], s2[k]) for k in s1), "sector gen not deterministic"
assert galaxy.generate_system(sid, star_idx, s1) == g2.generate_system(sid, star_idx, s2), \
    "system gen not deterministic"
assert galaxy.names == g2.names, "names not deterministic"
print("determinism: OK (fresh Galaxy instance reproduces identical output)")

# --- statistics across a radial sample of sectors -----------------------
sample = np.linspace(0, galaxy.n - 1, 12, dtype=int)
counts = collections.Counter()
n_hab = n_res = n_exo = 0
for s in sample:
    st = galaxy.generate_sector_stars(int(s))
    counts.update(CLASS_NAMES[c] for c in st["cls"])
    n_hab += st["habitable"].sum()
    n_res += st["resource_rich"].sum()
    n_exo += st["exotic"].sum()
tot = sum(counts.values())
print(f"\nclass mix over {len(sample)} sectors ({tot:,} stars):")
for k in CLASS_NAMES:
    print(f"  {k:>5}: {100.0 * counts.get(k, 0) / tot:6.3f} %")
print(f"stars per sector (avg over sample): {tot / len(sample):,.0f}")
print(f"flags per sector: habitable {n_hab / len(sample):.0f}, "
      f"resource {n_res / len(sample):.0f}, exotic {n_exo / len(sample):.0f}")

inner = galaxy.r < np.percentile(galaxy.r, 25)
outer = galaxy.r > np.percentile(galaxy.r, 75)
print(f"\ngradients: age inner {galaxy.age_gyr[inner].mean():.1f} / outer "
      f"{galaxy.age_gyr[outer].mean():.1f} Gyr | [Fe/H] inner "
      f"{galaxy.feh[inner].mean():+.2f} / outer {galaxy.feh[outer].mean():+.2f}")
