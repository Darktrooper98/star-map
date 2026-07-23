"""One-time tool: extract sector dot positions from the old RP galaxy map (IMG_7078.gif).

Each white dot on the old map is a sector (10,000 stars). This finds the dots,
filters out connection lines / nebula haze / text labels, and writes a catalog
to data/sectors_raw.json in normalized galaxy coordinates (core at origin,
map spans roughly [-150, 150] like the old star-map.py macro view).

Run with:  /Users/ashwinnimmal/opt/anaconda3/bin/python tools/extract_sectors.py
"""
import json
import os
import numpy as np
from PIL import Image
from scipy import ndimage

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GIF = os.path.join(ROOT, "IMG_7078.gif")
OUT = os.path.join(ROOT, "data", "sectors_raw.json")
MAP_HALF_EXTENT = 150.0  # match old star-map.py macro view scale


def main():
    img = Image.open(GIF).convert("RGB")
    a = np.asarray(img).astype(np.float32)
    h, w, _ = a.shape

    # Sector dots are bright and near-neutral (white), even inside tinted
    # faction territories. Lines are thin; nebulae are dim; text is white but
    # thin-stroked. Brightness threshold first:
    v = a.max(axis=2)
    sat = (a.max(axis=2) - a.min(axis=2)) / (a.max(axis=2) + 1e-6)

    def find_dots(mask):
        """Erode away thin structures (lines/text ~2-4 px; dots ~10+ px),
        then keep compact round blobs in a plausible dot-area band."""
        eroded = ndimage.binary_erosion(mask, iterations=3)
        labels, n = ndimage.label(eroded)
        if n == 0:
            return []
        idx = np.arange(1, n + 1)
        sizes = ndimage.sum(eroded, labels, index=idx)
        coms = ndimage.center_of_mass(eroded, labels, index=idx)
        slices = ndimage.find_objects(labels)
        out = []
        for i, (size, (cy, cx)) in enumerate(zip(sizes, coms)):
            sl = slices[i]
            bh = sl[0].stop - sl[0].start
            bw = sl[1].stop - sl[1].start
            if size < 8 or size > 400:
                continue
            aspect = max(bh, bw) / max(1, min(bh, bw))
            if aspect > 2.2:  # text fragments / line junctions are elongated
                continue
            fill = size / float(bh * bw)  # circles fill their bbox well
            if fill < 0.45:
                continue
            out.append((cx, cy, float(size)))
        return out

    # Pass 1: plain white dots (bright, near-neutral).
    dots = find_dots((v > 200) & (sat < 0.35))
    print(f"pass 1 (white dots): {len(dots)}")

    # Pass 2: dots tinted by faction territory overlays (bright but
    # saturated). Same shape filters keep out borders/nebulae; dedupe
    # against pass-1 dots.
    tinted = find_dots((v > 185) & (sat >= 0.35) & (sat < 0.85))
    from scipy.spatial import cKDTree
    if dots and tinted:
        tree = cKDTree([(d[0], d[1]) for d in dots])
        tinted = [t for t in tinted if not tree.query_ball_point((t[0], t[1]), r=9)]
    dots += tinted
    print(f"pass 2 (tinted dots): +{len(tinted)}, total {len(dots)}")

    # Map image coords -> galaxy coords: origin at image center (the core
    # black hole sits at the visual center), y up.
    cx0, cy0 = w / 2.0, h / 2.0
    sectors = []
    for cx, cy, size in dots:
        gx = (cx - cx0) / (w / 2.0) * MAP_HALF_EXTENT
        gy = -(cy - cy0) / (h / 2.0) * MAP_HALF_EXTENT
        sectors.append({"x": round(gx, 3), "y": round(gy, 3), "px_area": size})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump({"map_half_extent": MAP_HALF_EXTENT, "sectors": sectors}, f)
    print(f"wrote {len(sectors)} sectors -> {OUT}")

    # Diagnostic render for visual verification against the original.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs = [s["x"] for s in sectors]
    ys = [s["y"] for s in sectors]
    fig, ax = plt.subplots(figsize=(10, 10), facecolor="#050510")
    ax.set_facecolor("#050510")
    ax.scatter(xs, ys, s=6, c="white")
    ax.set_xlim(-155, 155); ax.set_ylim(-155, 155)
    ax.axis("off")
    fig.savefig(os.path.join(ROOT, "data", "extracted_preview.png"), dpi=110,
                facecolor="#050510", bbox_inches="tight")
    print("wrote data/extracted_preview.png")


if __name__ == "__main__":
    main()
