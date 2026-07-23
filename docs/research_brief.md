# Research Brief: Scaled-Down Stellar Forge for the RP Star Map

Consolidated E:D + astrophysics reference for implementation. Numbers are implementation targets, not
scholarly citations. Sources: klightspeed/EliteDangerousRegionMap (ID64 decode), 80.lv "Generating the
Universe in Elite: Dangerous" (Stellar Forge architecture), GDC talk recaps, EDSM statistics, standard
IMF/galactic-structure literature (Kroupa 2001, Chabrier 2003, Duchene & Kraus 2013, Kepler occurrence
papers, Milky Way disk structure surveys).

---

## TRACK A — How Elite: Dangerous / Stellar Forge actually works

### A1. Spatial hierarchy: sectors and boxels
- Galaxy volume is a fixed axis-aligned grid of **cubic sectors, 1280 ly per side**.
  Grid origin relative to Sol: `(x0, y0, z0) = (-49985, -40985, -24105)` ly (coords stored in 10-ly units).
- Each sector contains **8 overlapping octree layers** of "boxels", one per **mass code a–h**:

  | mass code | boxel side (ly) | boxels per sector | typical primary star |
  |---|---|---|---|
  | a | 10   | 128^3 = 2,097,152 | Y/T brown dwarfs |
  | b | 20   | 64^3              | L/T brown dwarfs, late M |
  | c | 40   | 32^3              | M, K dwarfs |
  | d | 80   | 16^3              | K, G, F, A dwarfs; white dwarfs; neutron stars |
  | e | 160  | 8^3               | A, F giants; B stars |
  | f | 320  | 4^3               | B stars, giants |
  | g | 640  | 2^3               | supergiants, small BH |
  | h | 1280 | 1                 | O, Wolf-Rayet, black holes, supergiants |

  Boxel side = `10 * 2^masscode` ly. Higher mass code = rarer, more massive systems; the octree layer
  IS the mass bucket. Every star system lives in exactly one boxel of one layer.

### A2. ID64 system address (the determinism backbone)
One 64-bit integer fully identifies a system; the seed IS the address. Bit layout (LSB first),
from klightspeed's RegionMap.py:

```
masscode  = id64 & 7                              # 3 bits: 0=a … 7=h
z_boxel   = (id64 >> 3)                & (0x3FFF >> mc)   # 14-mc bits
y_boxel   = (id64 >> (17 - mc))        & (0x1FFF >> mc)   # 13-mc bits
x_boxel   = (id64 >> (30 - 2*mc))      & (0x3FFF >> mc)   # 14-mc bits
n2        = remaining bits up to ~44             # serial index of system within boxel
bodyID    = top 9 bits (55-63)                   # body within system
coord_ly  = boxel_index << mc * 10 + origin      # e.g. x = ((xb << mc) * 10) - 49985
```
Design lesson: **address = seed = coordinates**. Given only the integer, the client regenerates the
whole system identically; nothing is stored. System names encode the same info
("Sector AB-C d3-42" = boxel letters + mass code + n2).

### A3. Galaxy-level inputs (what drives density/parameters per region)
- A **hand-authored 2D bitmap of the galaxy** (top-down star density image) extruded to 3D sets overall
  shape: spiral arms, bar/bulge. This is the only non-procedural galaxy input.
- Per-sector **"available mass" budget** + galaxy-scale **material (metallicity) and age distribution
  functions**: how much stuff is in a sector and how long it has been there.
- Generation is hierarchical: galaxy → sector/boxel → system → body; each child links to parent params.
- ~160,000 real systems (Hipparcos/Gliese/2MASS) are overlaid on top of procedural output (not needed for RP).

### A4. Star class / frequency in E:D
- Class chosen by sampling a mass from the boxel layer's mass bucket, then age + metallicity map
  mass → class (main sequence vs giant vs remnant): old + massive → remnant (WD/NS/BH), old + medium
  → giant, young + massive → O/B, etc.
- EDSM scanned-star mix (selection-biased but shows E:D targets, brown dwarfs included as primaries):
  M 36%, K 19%, L 10%, F 8%, G 6.5%, T 5%, Y 2.7%, T Tauri 2.5%, remainder A/B/O/giants/remnants.
  E:D roughly follows a real IMF but over-represents F relative to reality.

### A5. System contents from the seed
- E:D runs a coarse **accretion simulation**: total system mass + composition + angular momentum →
  single/binary/multiple star split; a protoplanetary disc is stepped through epochs; material
  aggregates into planets; effects modeled: solar wind stripping, tidal locking, gravitational heating,
  catastrophic events. Planet class = f(mass, composition, insolation).
- For the RP version a table-driven sampler (not a physics sim) reproduces the same outputs at ~1% of
  the cost — see recommendations.

---

## TRACK B — Astrophysics targets

### B1. Initial mass function (Kroupa 2001, piecewise power law dN/dm ∝ m^-α)
| mass range (M☉) | α |
|---|---|
| 0.01 – 0.08 (brown dwarfs) | 0.3 |
| 0.08 – 0.5 | 1.3 |
| 0.5 – 150 | 2.3 |

Sample by inverse-CDF over the three segments. (Chabrier alternative: lognormal m_c=0.079 M☉,
σ=0.69 below 1 M☉; power-law α=2.3 above — either is fine.)

### B2. Present-day stellar population (fraction of all stars, target table)
| type | fraction | notes |
|---|---|---|
| M V | 0.70 | red dwarfs dominate everywhere |
| K V | 0.11 | |
| G V | 0.06 | |
| F V | 0.030 | |
| A V | 0.006 | young-ish regions only (lifetime ~1 Gyr) |
| B V | 0.0012 | arms / young regions (lifetime <400 Myr) |
| O V | 0.00003 | only in active star-forming regions (lifetime <10 Myr) |
| White dwarf | 0.06 | old regions richer |
| Giants (K/M III) | 0.008 | old regions |
| Neutron star | 0.002 | |
| Black hole | 0.001 | |
| L/T/Y brown dwarfs | ~0.25 extra objects per star | include if E:D flavor wanted |

(Solar neighborhood observed: ~76% of main-sequence stars are M, ~12% K, ~7.6% G; WDs ~6% of all stars.)

### B3. Radial / vertical structure of a spiral galaxy
- Stellar density: `ρ(R,z) = ρ0 · exp(-R/h_R) · sech²(z/(2·h_z))` (or exp(-|z|/h_z)).
  Milky Way: h_R ≈ 2.4–3.9 kpc (adopt **2.6 kpc**); thin disk h_z ≈ 300 pc; thick disk h_z ≈ 900 pc
  with ~10% of stars, older/metal-poorer.
- Bulge/bar: central component, ~20–25% of stellar mass, radius ~2–3 kpc, roughly
  `ρ ∝ exp(-(R/1 kpc)²)` blob or short bar; stars old (>10 Gyr), high density (10–100× local).
- Spiral arms: **density perturbation on top of the disk**, not separate objects. Contrast vs interarm:
  ~1.5–3× for old stars, ~5–10× for gas/young stars. Model as logarithmic spirals
  `θ(R) = θ0 + ln(R/R0)/tan(p)`, pitch angle p ≈ 12–15°, 2 or 4 arms; multiply disk density by
  `1 + C·exp(-(Δθ/w)²)` with larger C for the young-star channel.
- Metallicity gradient: `[Fe/H](R) ≈ +0.3 − 0.06·(R − R_bulge) dex/kpc` (measured range −0.035 to
  −0.08 dex/kpc); scatter ±0.2 dex. Metallicity ↓ with height above plane and with age.
- Age gradient (inside-out formation): bulge/inner disk oldest (10–13 Gyr), solar radius mixed
  (mean ~5–7 Gyr), outer disk younger on average; BUT current star formation is concentrated in arms
  at all radii where gas exists (~3–15 kpc).

### B4. Where star formation happens (young stars, HII regions, nebulae)
- In arms within the gas-rich annulus (roughly 0.25–1.2 × R_solar-equivalent), plus scattered weak
  formation in interarm. Bulge and far outer disk: essentially none → no O/B, no T Tauri, no nebulae
  there; add extra WDs/NS/giants instead.
- Practical rule: give each sector a **star-formation-rate scalar** = arm boost × gas annulus ×
  (1 − bulge fraction); O/B/A/T-Tauri weights scale with it; remnant/giant weights scale with age.

### B5. Multiplicity and planets (gameplay flags)
| primary | binary/multiple fraction |
|---|---|
| O/B | 0.7–1.0 |
| A | 0.5–0.6 |
| F/G | 0.45 |
| K | 0.35 |
| M | 0.27 |
| brown dwarf | 0.15 |

- Planets: essentially all M/K/G/F stars host ≥1 planet; mean ~2.5 planets per M dwarf (P<200 d alone),
  similar or higher totals for FGK. Recommended: planet count ~ Poisson(λ), λ = 5–7 for single
  FGK/M, λ = 3 for close binaries, λ = 1 for remnants/OB.
- Habitable-zone terrestrial (η⊕): ~0.1–0.2 per M dwarf (conservative HZ), ~0.1–0.2 per FGK (large
  uncertainty). Recommended flag: P(HZ terrestrial) = 0.15 for F/G/K, 0.10 for M, 0.02 otherwise;
  boost ×1.5 at [Fe/H] > 0, ×0.5 at [Fe/H] < −0.5 (giant-planet occurrence rises steeply with
  metallicity; use a milder factor for terrestrials).
- HZ distance scales as `sqrt(L/L☉)` AU; gas giants more common around metal-rich and higher-mass stars.

---

## Recommended simplifications for a ~2000-sector RP galaxy

1. **Sector = the boxel.** Skip E:D's 8-layer octree; with only 2000 sectors × 10k stars, one layer
   suffices. Keep the *idea* of mass codes as a per-star "mass bucket" drawn from the IMF instead.
2. **Address = seed.** Pack `(sector_id, star_index)` into one 64-bit int
   (e.g. `sector_id << 32 | star_index`); hash it (e.g. SplitMix64) to seed a PRNG for on-demand,
   stateless, deterministic system generation. Never store generated systems.
3. **Galaxy shape from the GIF.** Load the old map image as the density bitmap (blur + threshold),
   exactly as Stellar Forge uses its 2D galaxy image; place the ~2000 sectors by importance-sampling
   or by reusing the white-dot positions directly. Add z by sampling `sech²(z/2h_z)` with h_z ≈ 4% of
   galaxy radius (thicker in bulge).
4. **Per-sector parameters (computed once, deterministically from sector position):**
   `R_frac` (0 at core, 1 at rim), `in_arm` (from bitmap brightness or log-spiral test),
   `age` = lerp(12 Gyr → 3 Gyr, core → rim) ± noise, `[Fe/H]` = +0.3 − 0.8·R_frac ± 0.2,
   `SFR` = arm × annulus factor (B4), `density` = bitmap value (display only).
5. **Star class table per sector, not per star physics.** Start from the B2 base table; reweight by
   sector age/SFR/metallicity (young: ×O/B/A/T-Tauri; old: ×WD/NS/BH/giants; then renormalize).
   Chi-by-eye targets beat running an IMF+stellar-evolution integration per star.
6. **System generation = table-driven, not accretion sim.** From the seeded PRNG: primary class →
   mass/temp/radius/luminosity from small per-class ranges; multiplicity roll (B5); planet count
   Poisson; each planet: type from (distance vs snow line, metallicity), simple orbit spacing
   (geometric ratio 1.4–2.0 from 0.05 AU), HZ flag via `sqrt(L)` scaling. This reproduces E:D-like
   output at trivial cost.
7. **Suggested global parameters:** galaxy radius = 1.0 (normalize); h_R = 0.25; bulge radius = 0.15
   with 20% of sectors; 2–4 arms, pitch 13°, arm density contrast 2× (old) / 6× (young weighting);
   metallicity gradient −0.8 dex core→rim; thin:thick disk 90:10 with h_z 0.03 : 0.09.
8. **Skip entirely:** real-star catalog overlay, per-planet accretion epochs, thick-disk kinematics,
   3D metallicity noise fields. None are visible at RP scale.
