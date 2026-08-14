"""Interactive RP galaxy map, powered by the scaled-down Stellar Forge
(stellarforge.py). Successor to star-map.py — same controls and features, but
the galaxy shape comes from the old hand-made map and everything is generated
deterministically: revisiting a sector or system always shows the same stars.

Views:
  Galaxy  — sectors (white dots), hyperlanes, nebulae, core black hole.
            WASD pan, scroll zoom, click a sector to enter it.
  Core    — the desolate galactic core around the TRF black hole.
  Sector  — that sector's stars (~80,000 on average) with class colors and
            the Habitable / Resource Rich / Exotic overlays. Click a star.
  System  — the star's full planetary system, generated on demand.

Run with:  /Users/ashwinnimmal/opt/anaconda3/bin/python starmap2.py
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import RadioButtons, Button
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
import numpy as np

from stellarforge import Galaxy, _rng, CLASS_NAMES

for _k in ('s', 'a', 'd', 'w'):
    for _map in (mpl.rcParams['keymap.save'], mpl.rcParams['keymap.pan']):
        if _k in _map:
            _map.remove(_k)

fig = plt.figure(figsize=(14, 8))
fig.patch.set_facecolor('#050510')
ui = {}
galaxy = Galaxy()


def add_scroll_zoom(ax, max_range=400.0, on_interact=None):
    def zoom(event):
        if event.inaxes != ax or event.xdata is None:
            return
        s = 1 / 1.15 if event.button == 'up' else 1.15 if event.button == 'down' else None
        if s is None:
            return
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        nw, nh = (x1 - x0) * s, (y1 - y0) * s
        if nw > max_range or nh > max_range:
            return
        rx = (x1 - event.xdata) / (x1 - x0)
        ry = (y1 - event.ydata) / (y1 - y0)
        if on_interact:
            on_interact()
        ax.set_xlim([event.xdata - nw * (1 - rx), event.xdata + nw * rx])
        ax.set_ylim([event.ydata - nh * (1 - ry), event.ydata + nh * ry])
        fig.canvas.draw_idle()
    ui['zoom_cid'] = fig.canvas.mpl_connect('scroll_event', zoom)


def add_wasd_pan(ax, on_interact=None):
    """Smooth held-key panning. Key press/release events maintain a set of
    held direction keys; a repeating timer glides the view a small step every
    tick while any key is down. This sidesteps the OS auto-repeat (initial
    delay + coarse repeats) that made single-press panning feel jumpy."""
    held = set()
    timer = fig.canvas.new_timer(interval=33)      # ~30 ticks/s
    step = 0.025                                    # view fraction per tick

    def tick():
        if not held:
            return
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        dx = (x1 - x0) * step * (('d' in held) - ('a' in held))
        dy = (y1 - y0) * step * (('w' in held) - ('s' in held))
        if on_interact:
            on_interact()
        ax.set_xlim(x0 + dx, x1 + dx)
        ax.set_ylim(y0 + dy, y1 + dy)
        fig.canvas.draw_idle()
    timer.add_callback(tick)
    ui['pan_timer'] = timer

    def on_press(event):
        k = (event.key or '').lower()
        if k in 'wasd' and k not in held:
            if not held:
                timer.start()
            held.add(k)

    def on_release(event):
        k = (event.key or '').lower()
        held.discard(k)
        if not held:
            timer.stop()

    ui['key_cid'] = fig.canvas.mpl_connect('key_press_event', on_press)
    ui['keyup_cid'] = fig.canvas.mpl_connect('key_release_event', on_release)


def defer_artist_during_interaction(artist, delay_ms=200, placeholder=None):
    """Returns an on_interact callback that hides `artist` while the view is
    being panned/zoomed and restores it shortly after interaction stops. If
    `placeholder` is given, it is shown while the heavy artist is hidden
    (e.g. a decimated star subset standing in for the full field).

    Rationale: expensive artists (viewport-sized image resampling, ~90k-marker
    scatters) make every pan/zoom frame slow. Swapping them out during
    interaction keeps frames fast; a single-shot timer brings them back.
    """
    timer = fig.canvas.new_timer(interval=delay_ms)
    timer.single_shot = True

    def restore():
        artist.set_visible(True)
        if placeholder is not None:
            placeholder.set_visible(False)
        fig.canvas.draw_idle()
    timer.add_callback(restore)
    ui['bg_timer'] = timer          # keep a reference so it isn't GC'd

    def on_interact():
        artist.set_visible(False)
        if placeholder is not None:
            placeholder.set_visible(True)
        timer.stop()
        timer.start()
    return on_interact


def clear_fig():
    for cid_key in ('pick_cid', 'zoom_cid', 'key_cid', 'keyup_cid'):
        if cid_key in ui:
            fig.canvas.mpl_disconnect(ui.pop(cid_key))
    for timer_key in ('bg_timer', 'pan_timer'):
        if timer_key in ui:
            ui[timer_key].stop()
    fig.clf()
    ui.clear()


def back_button(label, callback):
    ax_b = plt.axes([0.02, 0.85, 0.17, 0.05])
    btn = Button(ax_b, label, color='#2a2a3a', hovercolor='#3a3a4a')
    btn.label.set_color('white')
    btn.label.set_fontfamily('monospace')
    btn.label.set_weight('bold')
    btn.on_clicked(callback)
    ui['btn_back'] = btn


# --------------------------------------------------------------------- galaxy
_BG = {}


def _galaxy_background(extent=185.0, px=1000):
    """Rasterize the static decoration (dust haze + nebulae) into one RGBA
    image. These layers never change, but as live artists (~6k soft markers +
    ~380 circle patches) they dominated redraw time; as a single imshow they
    cost one draw call. Kept at modest resolution: resampling the image is
    itself per-frame work, and these are soft glows. Baked once and cached."""
    if 'img' not in _BG:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        f = Figure(figsize=(px / 100, px / 100), dpi=100)
        FigureCanvasAgg(f)
        bax = f.add_axes([0, 0, 1, 1])
        bax.set_xlim(-extent, extent)
        bax.set_ylim(-extent, extent)
        bax.axis('off')
        f.patch.set_alpha(0.0)
        bax.patch.set_alpha(0.0)

        x, y = galaxy.xy[:, 0], galaxy.xy[:, 1]
        haze = _rng("render", "haze")
        hx = np.repeat(x, 3) + haze.normal(0, 6, galaxy.n * 3)
        hy = np.repeat(y, 3) + haze.normal(0, 6, galaxy.n * 3)
        # Marker size (points^2) equivalent to a 10.4-data-unit diameter blob.
        pts_per_unit = px / (2 * extent) * 72.0 / 100.0
        bax.scatter(hx, hy, color='#191938', s=(10.4 * pts_per_unit) ** 2,
                    alpha=0.05, edgecolors='none')
        for sid in np.flatnonzero(galaxy.nebula):
            nrng = _rng("render", "nebula", int(sid))
            cx, cy = x[sid], y[sid]
            for _ in range(int(nrng.integers(5, 9))):
                qx, qy = cx + nrng.normal(0, 4), cy + nrng.normal(0, 4)
                rad = nrng.uniform(6, 14)
                bax.add_patch(Circle((qx, qy), rad, color='#330044', alpha=0.06, lw=0))
                bax.add_patch(Circle((qx, qy), rad * 0.7, color='#9900cc', alpha=0.04, lw=0))
            bax.add_patch(Circle((cx, cy), 5.0, color='#ff00ff', alpha=0.15, lw=0))
            bax.add_patch(Circle((cx, cy), 2.0, color='#e066ff', alpha=0.5))

        f.canvas.draw()
        _BG['img'] = np.asarray(f.canvas.buffer_rgba()).copy()
        _BG['extent'] = extent
    return _BG['img'], _BG['extent']


def show_galaxy_map(event=None):
    clear_fig()
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.85])
    ax.set_facecolor('#050510')
    ax.axis('off')
    ax.set_title('Macro Galactic View - WASD to Pan, Scroll to Zoom, Click a Sector',
                 color='white', fontsize=16, pad=20, fontfamily='monospace')
    ax.set_xlim(-150, 150)
    ax.set_ylim(-110, 110)
    ax.set_aspect('equal')

    x, y = galaxy.xy[:, 0], galaxy.xy[:, 1]

    # Static decoration (dust haze + nebulae), pre-baked to one image.
    img, ext = _galaxy_background()
    bg = ax.imshow(img, extent=(-ext, ext, -ext, ext), origin='upper',
                   interpolation='bilinear', zorder=0)
    on_interact = defer_artist_during_interaction(bg)

    segs = [[(x[i], y[i]), (x[j], y[j])] for i, j in galaxy.lanes]
    ax.add_collection(LineCollection(segs, colors='white', alpha=0.35,
                                     linewidths=0.5, zorder=1))
    dots = ax.scatter(x, y, color='white', s=14, picker=True, pickradius=4, zorder=2)

    # TRF: the core black hole.
    ax.scatter(0, 0, color='#ff7700', s=2500, alpha=0.15, zorder=7)
    ax.scatter(0, 0, color='#ffaa00', s=900, alpha=0.4, zorder=8)
    core = ax.scatter(0, 0, color='black', s=300, picker=True, zorder=9,
                      edgecolors='#ffcc00', linewidths=1.5)

    add_scroll_zoom(ax, 400.0, on_interact=on_interact)
    add_wasd_pan(ax, on_interact=on_interact)

    def on_pick(evt):
        if getattr(evt, 'mouseevent', None) and evt.mouseevent.button != 1:
            return
        if evt.artist == core:
            show_core_map()
        elif evt.artist == dots and len(evt.ind):
            show_sector_map(int(evt.ind[0]))
    ui['pick_cid'] = fig.canvas.mpl_connect('pick_event', on_pick)
    fig.canvas.draw_idle()


# ----------------------------------------------------------------------- core
def show_core_map(event=None):
    clear_fig()
    rng = _rng("core", "stars")
    n = 900
    x, y = rng.random(n) * 100, rng.random(n) * 100
    keep = np.hypot(x - 50, y - 50) > 25
    x, y = x[keep], y[keep]

    ax = fig.add_axes([0.22, 0.05, 0.75, 0.85])
    ax.set_facecolor('#050510')
    ax.axis('off')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_title(f'Desolate Galactic Core - TRF ({len(x):,} ancient stars)',
                 color='white', fontsize=16, fontfamily='monospace')
    # The core hosts only ancient red stars and remnants.
    colors = np.array(['#ff4d4d', '#ffb05b', '#cfe0ff'])[rng.choice(3, len(x), p=[0.6, 0.25, 0.15])]
    ax.scatter(x, y, s=rng.random(len(x)) * 2 + 0.3, c=colors, alpha=0.7, edgecolors='none')
    ax.scatter(50, 50, s=40000, color='#ff7700', alpha=0.08, zorder=8)
    ax.scatter(50, 50, s=15000, color='#ffaa00', alpha=0.15, zorder=9)
    ax.scatter(50, 50, s=5000, color='black', edgecolors='white', linewidths=2, zorder=10)
    add_scroll_zoom(ax, 300.0)
    back_button('< Galaxy', show_galaxy_map)
    fig.canvas.draw_idle()


# --------------------------------------------------------------------- sector
def show_sector_map(sector_id, event=None):
    clear_fig()
    stars = galaxy.generate_sector_stars(sector_id)
    n = len(stars['x'])

    ax = fig.add_axes([0.22, 0.05, 0.75, 0.85])
    ax.set_facecolor('#050510')
    ax.axis('off')
    ax.set_xlim(-3, 103)
    ax.set_ylim(-3, 103)

    info = (f"{galaxy.names[sector_id]}  |  {n:,} stars  |  "
            f"age {galaxy.age_gyr[sector_id]:.1f} Gyr  |  "
            f"[Fe/H] {galaxy.feh[sector_id]:+.2f}")
    title = ax.set_title(info + '  (Click a star)', color='white',
                         fontsize=13, fontfamily='monospace')

    scatter = ax.scatter(stars['x'], stars['y'], s=stars['size'], c=stars['color'],
                         alpha=0.6, edgecolors='none', picker=True, pickradius=3)
    # Decimated stand-in shown while zooming (full field restores on idle).
    DEC = 6
    preview = ax.scatter(stars['x'][::DEC], stars['y'][::DEC],
                         s=stars['size'][::DEC] * 1.5, c=stars['color'][::DEC],
                         alpha=0.7, edgecolors='none', visible=False)
    # Filter modes swap the full field out for this scatter, which holds ONLY
    # the matching stars - filtered-out stars are gone, not dimmed, so they
    # can't be drawn over or picked instead of the dot the user aimed at.
    filt = ax.scatter([], [], visible=False, edgecolors='none',
                      picker=True, pickradius=5)
    state = {'mode': 'Star Class', 'filt_idx': np.empty(0, int)}

    full_interact = defer_artist_during_interaction(scatter, placeholder=preview)
    add_scroll_zoom(ax, 300.0,
                    on_interact=lambda: state['mode'] == 'Star Class' and full_interact())

    ax_radio = plt.axes([0.02, 0.40, 0.17, 0.25], facecolor='#0a0a1a')
    radio = RadioButtons(ax_radio, ('Star Class', 'Habitable Systems',
                                    'Resource Rich', 'Exotic Stars'),
                         activecolor='#00ffcc')
    for lab in radio.labels:
        lab.set_color('white')
        lab.set_fontsize(10)
        lab.set_fontfamily('monospace')

    overlays = {'Habitable Systems': ('habitable', '#00ffcc', 14.0),
                'Resource Rich': ('resource_rich', '#ffd700', 8.0),
                'Exotic Stars': ('exotic', '#ff00ff', 20.0)}

    def set_mode(label):
        state['mode'] = label
        if label in overlays:
            key, col, sz = overlays[label]
            idx = np.flatnonzero(stars[key])
            state['filt_idx'] = idx
            filt.set_offsets(np.column_stack((stars['x'][idx], stars['y'][idx])))
            filt.set_facecolors(col)
            filt.set_sizes(np.full(idx.size, sz))
            scatter.set_visible(False)
            preview.set_visible(False)
            filt.set_visible(True)
            title.set_text(f"{galaxy.names[sector_id]}  |  {label}: "
                           f"{idx.size:,} of {n:,} stars  (Click a star)")
        else:
            filt.set_visible(False)
            preview.set_visible(False)
            scatter.set_visible(True)
            title.set_text(info + '  (Click a star)')
        fig.canvas.draw_idle()
    radio.on_clicked(set_mode)
    ui['radio'] = radio

    def nearest_to_click(evt, xs, ys):
        """Of all candidate indices in the pick event, take the star closest
        to the actual click point (evt.ind[0] is just the lowest index and is
        usually the wrong star in an 80k-star field)."""
        ii = np.asarray(evt.ind)
        mx, my = evt.mouseevent.xdata, evt.mouseevent.ydata
        return int(ii[np.argmin((xs[ii] - mx) ** 2 + (ys[ii] - my) ** 2)])

    def on_pick(evt):
        if not len(getattr(evt, 'ind', [])) or evt.mouseevent.button != 1 \
                or evt.mouseevent.xdata is None:
            return
        if evt.artist == scatter and scatter.get_visible():
            show_system_view(sector_id, nearest_to_click(evt, stars['x'], stars['y']),
                             stars)
        elif evt.artist == filt and filt.get_visible():
            idx = state['filt_idx']
            sub = nearest_to_click(evt, stars['x'][idx], stars['y'][idx])
            show_system_view(sector_id, int(idx[sub]), stars)
    ui['pick_cid'] = fig.canvas.mpl_connect('pick_event', on_pick)

    back_button('< Galaxy', show_galaxy_map)
    fig.canvas.draw_idle()


# --------------------------------------------------------------------- system
def show_system_view(sector_id, star_index, stars):
    sysd = galaxy.generate_system(sector_id, star_index, stars)
    clear_fig()

    ax = fig.add_axes([0.30, 0.05, 0.65, 0.85])
    ax.set_facecolor('#050510')
    ax.axis('off')
    ax.set_aspect('equal')
    ax.set_title(f"{sysd['name']}  -  class {sysd['class']}  "
                 f"({sysd['mass_msun']} Msun, {sysd['age_gyr']} Gyr)",
                 color='white', fontsize=13, fontfamily='monospace')

    cls_idx = CLASS_NAMES.index(sysd['class'])
    star_col = ['#5b82ff', '#9bb0ff', '#ffffff', '#fff4ea', '#ffe56f', '#ffb05b',
                '#ff4d4d', '#cfe0ff', '#ffcc88', '#b0fff6', '#111118'][cls_idx]

    a_list = [p['a_au'] for p in sysd['planets']] or [1.0]
    lim = max(a_list) * 1.25
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    ax.scatter(0, 0, s=500, color=star_col, zorder=5,
               edgecolors='#ffcc00' if sysd['class'] == 'BH' else 'none',
               linewidths=1.5)
    if sysd['hz_au'] > 0.01 and sysd['class'] not in ('NS', 'BH', 'WD'):
        ax.add_patch(Circle((0, 0), sysd['hz_au'], fill=False, color='#00ffcc',
                            alpha=0.25, ls='--', lw=1.0))

    prng = _rng("render", "system", sector_id, star_index)
    type_cols = {'scorched rocky': '#c08050', 'terrestrial': '#88bb77',
                 'terrestrial (habitable)': '#00ffcc', 'gas dwarf': '#ddaa99',
                 'gas giant': '#ddbb66', 'ice giant': '#88bbee', 'icy body': '#cceeff'}
    for p in sysd['planets']:
        ax.add_patch(Circle((0, 0), p['a_au'], fill=False, color='white',
                            alpha=0.15, lw=0.6))
        ang = prng.uniform(0, 2 * np.pi)
        px, py = p['a_au'] * np.cos(ang), p['a_au'] * np.sin(ang)
        ax.scatter(px, py, s=np.clip(p['mass_earth'], 1, 400) ** 0.5 * 8 + 10,
                   color=type_cols.get(p['type'], 'white'), zorder=6)

    lines = [f"luminosity {sysd['luminosity_lsun']} Lsun",
             f"habitable zone ~{sysd['hz_au']} AU"]
    for c in sysd['companions']:
        lines.append(f"companion star: {c['mass']} Msun @ {c['sep_au']} AU")
    lines.append("")
    for p in sysd['planets']:
        lines.append(f"{p['name'].split()[-1]}: {p['type']}")
        lines.append(f"   {p['a_au']} AU, {p['mass_earth']} Me")
    if not sysd['planets']:
        lines.append("no planets")
    fig.text(0.02, 0.72, "\n".join(lines), color='#aaccdd', fontsize=9,
             fontfamily='monospace', va='top')

    add_scroll_zoom(ax, lim * 4)
    back_button('< Sector', lambda e: show_sector_map(sector_id))
    fig.canvas.draw_idle()


if __name__ == '__main__':
    show_galaxy_map()
    plt.show()
