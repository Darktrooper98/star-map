"""Scaled-down Stellar Forge: deterministic procedural galaxy for the RP star map.

Hierarchy (mirroring Elite: Dangerous' approach, scaled down):
  Galaxy -> Sectors (from the old hand-made map; each = 10,000 stars) -> Boxels
  (spatial sub-cells of a sector) -> Stars -> Systems (planets, on demand).

Every level is generated deterministically from GALAXY_SEED, so any sector or
individual star system can be regenerated on demand and is always identical.

Each sector carries astrophysical parameters derived from its position in the
galaxy (distance from core, local arm density):
  age_gyr          mean stellar age (inside-out formation: old core, young rim)
  feh              metallicity [Fe/H] (radial gradient, ~-0.9 dex core->rim)
  young_frac       fraction of recently formed stars (high in dense arm regions)
  nebula           True for the most actively star-forming sectors

Star generation samples masses from the Kroupa IMF, evolves them against the
star's age (main-sequence lifetime ~ 10*M^-2.5 Gyr) into giants / white dwarfs /
neutron stars / black holes, and derives gameplay flags (habitable, resource
rich, exotic). generate_system() expands any star into a full planetary system.
"""
import hashlib
import json
import os

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import coo_matrix

GALAXY_SEED = "TRF-GALAXY-V1"
ROOT = os.path.dirname(os.path.abspath(__file__))
SECTOR_FILE = os.path.join(ROOT, "data", "sectors_raw.json")
STARS_PER_SECTOR = 10_000
NUM_NEBULAE = 25
BOXELS_PER_SIDE = 8          # sector view is an 8x8 boxel grid

# Spectral classes, coldest last. Mass ranges in solar masses (main sequence).
CLASS_NAMES = ["O", "B", "A", "F", "G", "K", "M", "WD", "GIANT", "NS", "BH"]
CLASS_COLORS = ["#5b82ff", "#9bb0ff", "#ffffff", "#fff4ea", "#ffe56f",
                "#ffb05b", "#ff4d4d", "#cfe0ff", "#ffcc88", "#b0fff6", "#3a3a4a"]
MS_MASS_EDGES = [16.0, 2.1, 1.4, 1.04, 0.8, 0.45, 0.08]  # lower edge of O..M

_SYL_A = ["Au", "Bei", "Cro", "Dra", "Eol", "Fel", "Gru", "Hyo", "Ithe", "Jen",
          "Kai", "Lyra", "Mor", "Nyx", "Oph", "Pra", "Qui", "Rhe", "Sag", "Tau",
          "Ur", "Vel", "Wre", "Xi", "Yso", "Zel"]
_SYL_B = ["ba", "ce", "di", "fa", "ga", "he", "ki", "la", "me", "ni",
          "pa", "ra", "se", "ta", "ve", "wa", "xe", "ya", "zo", "ru"]
_SYL_C = ["", "", "n", "s", "th", "x", "m", "r", "l", "k"]


def _seed(*parts):
    """Stable 64-bit seed from a path of identifiers (Python's hash() is salted)."""
    key = "/".join(str(p) for p in (GALAXY_SEED,) + parts)
    return int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")


def _rng(*parts):
    return np.random.default_rng(_seed(*parts))


def _sector_name(rng):
    n = _SYL_A[rng.integers(len(_SYL_A))] + _SYL_B[rng.integers(len(_SYL_B))]
    if rng.random() < 0.5:
        n += _SYL_B[rng.integers(len(_SYL_B))]
    n += _SYL_C[rng.integers(len(_SYL_C))]
    return n


def _kroupa_masses(rng, n, m_max=100.0):
    """Sample n masses [0.08, m_max] Msun from the Kroupa (2001) broken power
    law: dN/dM ~ M^-1.3 (0.08-0.5), M^-2.3 (>0.5)."""
    def seg_integral(a, lo, hi):
        return (hi ** (1 - a) - lo ** (1 - a)) / (1 - a)

    w1 = seg_integral(1.3, 0.08, 0.5)
    w2 = 0.5 ** (2.3 - 1.3) * seg_integral(2.3, 0.5, m_max)  # continuity at 0.5
    u = rng.random(n)
    pick_low = u < w1 / (w1 + w2)
    m = np.empty(n)

    def inv_cdf(a, lo, hi, size, r):
        t = r.random(size)
        return (lo ** (1 - a) + t * (hi ** (1 - a) - lo ** (1 - a))) ** (1 / (1 - a))

    m[pick_low] = inv_cdf(1.3, 0.08, 0.5, pick_low.sum(), rng)
    m[~pick_low] = inv_cdf(2.3, 0.5, m_max, (~pick_low).sum(), rng)
    return m


class Galaxy:
    """Sector catalog with derived astrophysical parameters and lane network."""

    def __init__(self, sector_file=SECTOR_FILE):
        with open(sector_file) as f:
            raw = json.load(f)
        pts = np.array([[s["x"], s["y"]] for s in raw["sectors"]])
        # Stable ordering (and therefore stable sector ids/names) regardless of
        # extraction order: sort by angle then radius.
        r = np.hypot(pts[:, 0], pts[:, 1])
        order = np.lexsort((r, np.arctan2(pts[:, 1], pts[:, 0])))
        pts = pts[order]
        self.xy = pts
        self.n = len(pts)
        self.half_extent = raw["map_half_extent"]

        self.r = np.hypot(pts[:, 0], pts[:, 1])
        r_norm = self.r / self.r.max()

        # Local density from mean distance to 6 nearest neighbours.
        tree = cKDTree(pts)
        d, _ = tree.query(pts, k=7)
        mean_nn = d[:, 1:].mean(axis=1)
        dens = 1.0 / np.maximum(mean_nn, 1e-3) ** 2
        self.density = (dens - dens.min()) / (dens.max() - dens.min())

        # Per-sector parameters, deterministic noise per sector id.
        noise = np.array([_rng("sector", i, "params").normal(0, 1, 3)
                          for i in range(self.n)])
        self.age_gyr = np.clip(11.5 - 6.5 * r_norm + noise[:, 0] * 1.2, 0.3, 13.0)
        self.feh = np.clip(0.35 - 0.9 * r_norm + noise[:, 1] * 0.12, -1.6, 0.6)
        self.young_frac = np.clip(0.03 + 0.30 * self.density
                                  + 0.10 * r_norm + noise[:, 2] * 0.05, 0.0, 0.55)

        # Nebulae: highest-SFR sectors, greedily thinned so they scatter
        # across the arms instead of chaining along one dense region.
        sfr = self.young_frac * (0.5 + self.density)
        self.nebula = np.zeros(self.n, bool)
        min_sep = self.r.max() * 0.12
        for i in np.argsort(sfr)[::-1]:
            if self.nebula.sum() >= NUM_NEBULAE:
                break
            picked = np.flatnonzero(self.nebula)
            if picked.size == 0 or np.min(np.hypot(
                    *(self.xy[picked] - self.xy[i]).T)) > min_sep:
                self.nebula[i] = True

        self.names = self._unique_names()
        self.lanes = self._build_lanes(tree)

    def _unique_names(self):
        names, used = [], set()
        for i in range(self.n):
            rng = _rng("sector", i, "name")
            name = _sector_name(rng)
            while name in used:
                name = _sector_name(rng)
            used.add(name)
            code = f"{chr(65 + int(rng.integers(26)))}{chr(65 + int(rng.integers(26)))}"
            names.append(f"{name} {code}-{i % 10}")
        return names

    def _build_lanes(self, tree):
        """Hyperlane network: MST over the k-NN graph (guarantees connectivity)
        plus short local links, like the old map's web of lanes."""
        k = min(6, self.n - 1)
        d, idx = tree.query(self.xy, k=k + 1)
        rows = np.repeat(np.arange(self.n), k)
        cols = idx[:, 1:].ravel()
        vals = d[:, 1:].ravel()
        mst = minimum_spanning_tree(
            coo_matrix((vals, (rows, cols)), shape=(self.n, self.n)))
        pairs = {tuple(sorted(p)) for p in zip(*mst.nonzero())}
        # Extra local links for the web-like look, but only short ones.
        thresh = np.median(d[:, 1]) * 2.2
        for i in range(self.n):
            for dist, j in zip(d[i, 1:4], idx[i, 1:4]):
                if dist < thresh:
                    pairs.add(tuple(sorted((i, int(j)))))
        return sorted(pairs)

    # ------------------------------------------------------------------ stars
    def generate_sector_stars(self, sector_id, n_stars=STARS_PER_SECTOR):
        """Deterministically generate the stars of one sector.

        Returns a dict of arrays (all length n_stars): x, y (0-100 sector
        coords), cls (index into CLASS_NAMES), mass, age_gyr, color, size,
        habitable, resource_rich, exotic.
        """
        rng = _rng("sector", sector_id, "stars")

        # --- positions: boxel grid with lognormal weights (Stellar Forge's
        # boxels, scaled down) + per-boxel jitter -> clumpy, stable texture.
        nb = BOXELS_PER_SIDE
        w = _rng("sector", sector_id, "boxels").lognormal(0.0, 0.8, nb * nb)
        w /= w.sum()
        boxel = rng.choice(nb * nb, size=n_stars, p=w)
        bx, by = boxel % nb, boxel // nb
        side = 100.0 / nb
        # Soft boxel edges: gaussian scatter past the cell bounds so the grid
        # doesn't read as blocks.
        x = np.clip((bx + rng.random(n_stars)) * side + rng.normal(0, side * 0.25, n_stars), 0, 100)
        y = np.clip((by + rng.random(n_stars)) * side + rng.normal(0, side * 0.25, n_stars), 0, 100)

        # --- ages: mix of the sector's old population and recent formation.
        young = rng.random(n_stars) < self.young_frac[sector_id]
        age = np.where(young,
                       rng.random(n_stars) * 2.0,
                       np.clip(rng.normal(self.age_gyr[sector_id], 1.0, n_stars),
                               0.2, 13.2))

        # --- masses from the IMF; evolve against age.
        mass = _kroupa_masses(rng, n_stars)
        t_ms = 10.0 * mass ** -2.5                      # MS lifetime, Gyr
        cls = np.searchsorted(-np.array(MS_MASS_EDGES), -mass)  # 0=O .. 6=M
        cls = np.clip(cls, 0, 6)

        evolved = age > t_ms
        giant = (~evolved) & (age > 0.85 * t_ms)
        cls[giant] = 8                                   # GIANT
        cls[evolved & (mass < 8)] = 7                    # WD
        cls[evolved & (mass >= 8) & (mass < 25)] = 9     # NS
        cls[evolved & (mass >= 25)] = 10                 # BH

        # --- gameplay flags. Habitability needs a stable MS F/G/K (or lucky M)
        # host and is boosted by metallicity (planet occurrence correlation).
        z_boost = np.clip(10 ** (0.6 * self.feh[sector_id]), 0.15, 2.5)
        p_hab = np.select([np.isin(cls, [3, 4, 5]), cls == 6],
                          [0.02, 0.004], default=0.0) * z_boost
        habitable = rng.random(n_stars) < p_hab
        resource_rich = rng.random(n_stars) < 0.03 * np.clip(z_boost, 0.3, 2.0)
        exotic = np.isin(cls, [9, 10]) | (rng.random(n_stars) < 0.001)

        color = np.array(CLASS_COLORS)[cls]
        size = np.clip(mass, 0.08, 30) ** 0.5 * 2.4
        size[np.isin(cls, [7, 9])] = 1.0                 # compact remnants
        size[cls == 10] = 2.5

        return {"x": x, "y": y, "cls": cls, "mass": mass, "age_gyr": age,
                "color": color, "size": size, "habitable": habitable,
                "resource_rich": resource_rich, "exotic": exotic}

    # ---------------------------------------------------------------- systems
    def generate_system(self, sector_id, star_index, stars=None):
        """Expand one star into a full system, deterministically.

        `stars` may pass a cached generate_sector_stars() result to avoid
        regenerating the sector.
        """
        if stars is None:
            stars = self.generate_sector_stars(sector_id)
        rng = _rng("sector", sector_id, "system", star_index)
        cls = int(stars["cls"][star_index])
        mass = float(stars["mass"][star_index])
        flagged_hab = bool(stars["habitable"][star_index])

        name = f"{self.names[sector_id]} {star_index}"
        lum = np.clip(mass, 0.08, 60) ** 3.5             # L/Lsun, crude MS law
        if cls == 7:
            lum = 0.001
        elif cls == 8:
            lum = lum * 40
        elif cls in (9, 10):
            lum = 0.0001

        companions = []
        p_multi = [0.8, 0.7, 0.6, 0.5, 0.45, 0.3, 0.25, 0.1, 0.4, 0.1, 0.2][cls]
        if rng.random() < p_multi:
            companions.append({"mass": round(mass * float(rng.uniform(0.15, 0.9)), 3),
                               "sep_au": round(float(rng.lognormal(3.0, 1.5)), 2)})

        # Planet count rises with metallicity; dead remnants keep few.
        lam = 4.0 * np.clip(10 ** (0.5 * self.feh[sector_id]), 0.2, 2.5)
        if cls >= 7:
            lam *= 0.4
        n_planets = int(rng.poisson(lam))
        if flagged_hab:
            n_planets = max(n_planets, 1)

        hz_center = float(np.sqrt(max(lum, 1e-4)))       # ~1 AU * sqrt(L)
        planets = []
        a = 0.05 * max(hz_center, 0.2)
        hab_slot = rng.integers(n_planets) if flagged_hab else -1
        for p in range(n_planets):
            a *= float(rng.uniform(1.4, 2.1))            # spacing law
            if p == hab_slot:
                a = float(hz_center * rng.uniform(0.95, 1.25))
            rel = a / max(hz_center, 1e-3)
            if p == hab_slot:
                ptype, pm = "terrestrial (habitable)", float(rng.uniform(0.5, 2.0))
            elif rel < 0.5:
                ptype, pm = "scorched rocky", float(rng.uniform(0.05, 1.5))
            elif rel < 1.6:
                ptype = "terrestrial" if rng.random() < 0.6 else "gas dwarf"
                pm = float(rng.uniform(0.1, 4.0))
            elif rel < 6.0:
                # Gas giants are rare around low-mass hosts.
                if rng.random() < min(1.0, mass / 0.5):
                    ptype, pm = "gas giant", float(rng.uniform(20, 600))
                else:
                    ptype, pm = "icy body", float(rng.uniform(0.05, 5.0))
            else:
                ptype, pm = "ice giant" if rng.random() < 0.6 else "icy body", \
                            float(rng.uniform(0.05, 40))
            planets.append({"name": f"{name} {chr(98 + p)}", "a_au": round(a, 3),
                            "type": ptype, "mass_earth": round(pm, 2)})

        return {"name": name, "sector": self.names[sector_id],
                "class": CLASS_NAMES[cls], "mass_msun": round(mass, 3),
                "age_gyr": round(float(stars["age_gyr"][star_index]), 2),
                "luminosity_lsun": round(float(lum), 4),
                "hz_au": round(hz_center, 3), "companions": companions,
                "planets": planets}


if __name__ == "__main__":
    g = Galaxy()
    print(f"{g.n} sectors, {len(g.lanes)} lanes, {g.nebula.sum()} nebulae")
    s = g.generate_sector_stars(100)
    import collections
    counts = collections.Counter(CLASS_NAMES[c] for c in s["cls"])
    print("sector 100:", g.names[100], "| class mix:", dict(counts))
    print("habitable:", int(s['habitable'].sum()),
          "resource:", int(s['resource_rich'].sum()),
          "exotic:", int(s['exotic'].sum()))
    sys0 = g.generate_system(100, int(np.flatnonzero(s["habitable"])[0]), s)
    print(json.dumps(sys0, indent=1))
