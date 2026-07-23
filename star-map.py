import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import RadioButtons, Button
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
import numpy as np
from scipy.spatial import KDTree
import heapq

# Disable Matplotlib's default 's' key shortcut for saving figures so WASD works cleanly
if 's' in mpl.rcParams['keymap.save']:
    mpl.rcParams['keymap.save'].remove('s')

# Set up the main application window
fig = plt.figure(figsize=(14, 8))
fig.patch.set_facecolor('#050510')

# Dictionary to hold UI elements
ui_elements = {}

def add_scroll_zoom(ax, max_range=400.0):
    """Enables smooth scroll zooming centered on the mouse with a maximum zoom-out limit."""
    def zoom_fun(event):
        if event.inaxes != ax:
            return
        
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        
        base_scale = 1.15
        if event.button == 'up':
            scale_factor = 1.0 / base_scale
        elif event.button == 'down':
            scale_factor = base_scale
        else:
            return

        xdata = event.xdata
        ydata = event.ydata
        
        if xdata is None or ydata is None:
            return

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
        
        if new_width > max_range or new_height > max_range:
            return

        relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])
        
        ax.set_xlim([xdata - new_width * (1 - relx), xdata + new_width * relx])
        ax.set_ylim([ydata - new_height * (1 - rely), ydata + new_height * rely])
        ax.figure.canvas.draw_idle()

    fig.canvas.mpl_connect('scroll_event', zoom_fun)

def show_galaxy_map(event=None):
    """Draws the macro map with organic multi-layered nebula formations (capped at 25)."""
    fig.clf() 
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.85])
    ax.set_facecolor('#050510')
    ax.axis('off')
    ax.set_title('Macro Galactic View - WASD to Pan, Scroll to Zoom', 
                 color='white', fontsize=18, pad=20, fontfamily='monospace')
    
    ax.set_xlim(-150, 150)
    ax.set_ylim(-110, 110)

    # --- 1. Generate Nebula and Dust Structure ---
    neb_theta = np.linspace(0, 2 * np.pi * 3, 300) 
    neb_r = (neb_theta ** 1.2 * 12) + 12  
    neb_offset = (np.arange(300) % 7) * (2 * np.pi / 7)
    
    scale_neb = neb_r * 0.1
    neb_x = neb_r * np.cos(neb_theta + neb_offset) + np.array([np.random.normal(0, s) for s in scale_neb])
    neb_y = neb_r * np.sin(neb_theta + neb_offset) + np.array([np.random.normal(0, s) for s in scale_neb])
    ax.scatter(neb_x, neb_y, color='#110033', s=12000, alpha=0.08, zorder=0)

    # --- 2. Generate Star Nodes with Tighter Core Clearance ---
    num_nodes = np.random.randint(1800, 2200)
    arms = 5
    theta = np.linspace(0, 2 * np.pi * 2, num_nodes) 
    
    r = (theta ** 1.3 * 10) + 25 
    arm_offset = (np.arange(num_nodes) % arms) * (2 * np.pi / arms)
    
    scale_nodes = r * 0.08
    x = r * np.cos(theta + arm_offset) + np.array([np.random.normal(0, s) for s in scale_nodes])
    y = r * np.sin(theta + arm_offset) + np.array([np.random.normal(0, s) for s in scale_nodes])
    
    core_distances = np.sqrt(x**2 + y**2)
    min_core_clearance = 15.0 
    valid_macro_nodes = core_distances > min_core_clearance
    x = x[valid_macro_nodes]
    y = y[valid_macro_nodes]
    num_nodes = len(x)

    nodes_scatter = ax.scatter(x, y, color='white', s=15, alpha=1.0, picker=True, zorder=2)
    
    # --- 3. Core-Connected Network Generation with Outlier Sweep ---
    points = np.column_stack((x, y))
    core_anchor_idx = np.argmin(x**2 + y**2)
    
    tree = KDTree(points)
    pairs = set()
    
    core_connections = 0
    max_core_lanes = 5
    
    visited = [False] * num_nodes
    min_heap = []
    
    visited[core_anchor_idx] = True
    dists, indices = tree.query(points[core_anchor_idx], k=min(15, num_nodes))
    for d, idx in zip(dists[1:], indices[1:]):
        heapq.heappush(min_heap, (d, core_anchor_idx, idx))
        
    while min_heap:
        d, u, v = heapq.heappop(min_heap)
        
        if visited[v]:
            if d < 30 and len(pairs) < num_nodes * 3:
                pair = tuple(sorted((u, v)))
                if core_anchor_idx in pair:
                    if core_connections < max_core_lanes and pair not in pairs:
                        pairs.add(pair)
                        core_connections += 1
                else:
                    pairs.add(pair)
            continue
            
        visited[v] = True
        
        pair = tuple(sorted((u, v)))
        if core_anchor_idx in pair:
            if core_connections < max_core_lanes:
                pairs.add(pair)
                core_connections += 1
        else:
            pairs.add(pair)
            
        dists_v, indices_v = tree.query(points[v], k=min(7, num_nodes))
        for dv, iv in zip(dists_v[1:], indices_v[1:]):
            if not visited[iv]:
                heapq.heappush(min_heap, (dv, v, iv))

    for i in range(num_nodes):
        dists, indices = tree.query(points[i], k=5)
        for d, j in zip(dists[1:], indices[1:]):
            if d < 25:
                pair = tuple(sorted((i, j)))
                if core_anchor_idx in pair:
                    if core_connections < max_core_lanes:
                        pairs.add(pair)
                        core_connections += 1
                else:
                    pairs.add(pair)

    connected_nodes = set(node for p in pairs for node in p)
    for i in range(num_nodes):
        if i not in connected_nodes:
            dists, indices = tree.query(points[i], k=2)
            nearest = indices[1]
            pair = tuple(sorted((i, nearest)))
            if core_anchor_idx in pair:
                if core_connections < max_core_lanes:
                    pairs.add(pair)
                    core_connections += 1
            else:
                pairs.add(pair)

    segments = [[(x[i], y[i]), (x[j], y[j])] for i, j in pairs]
    lc = LineCollection(segments, colors='white', alpha=0.5, linewidths=0.6, zorder=1)
    ax.add_collection(lc)

    # --- 4. Multi-Layered Organic Nebula Sectors (Capped at 25) ---
    num_clusters = min(25, num_nodes)
    cluster_indices = np.random.choice(num_nodes, size=num_clusters, replace=False)
    
    np.random.seed(42) # Keep layout stable across redraws if desired, or remove for random
    for idx in cluster_indices:
        cx, cy = x[idx], y[idx]
        
        # Draw overlapping wisps of gas to create a real nebula cloud look
        num_puffs = np.random.randint(5, 9)
        for _ in range(num_puffs):
            puff_x = cx + np.random.normal(0, 4.0)
            puff_y = cy + np.random.normal(0, 4.0)
            puff_radius = np.random.uniform(6.0, 14.0)
            
            # Soft outer gas layer
            ax.add_patch(Circle((puff_x, puff_y), radius=puff_radius, 
                                color='#330044', alpha=0.06, linewidth=0, zorder=3))
            # Mid layer magenta glow
            ax.add_patch(Circle((puff_x, puff_y), radius=puff_radius * 0.7, 
                                color='#9900cc', alpha=0.04, linewidth=0, zorder=4))

        # Core bright glowing center of the nebula cloud
        ax.add_patch(Circle((cx, cy), radius=5.0, color='#ff00ff', alpha=0.15, linewidth=0, zorder=5))
        ax.add_patch(Circle((cx, cy), radius=2.0, color='#e066ff', alpha=0.5, 
                            edgecolor='white', linewidth=0.6, zorder=6))
    
    # --- 5. Core Black Hole and Void ---
    ax.scatter(0, 0, color='#ff7700', s=2500, alpha=0.15, zorder=7) 
    ax.scatter(0, 0, color='#ffaa00', s=900, alpha=0.4, zorder=8)   
    core_node = ax.scatter(0, 0, color='black', s=300, picker=True, zorder=9, edgecolors='#ffcc00', linewidths=1.5)
    
    add_scroll_zoom(ax, max_range=400.0)

    # --- 6. WASD Panning Support ---
    def on_key(event):
        if event.key not in ['w', 'a', 's', 'd', 'W', 'A', 'S', 'D']:
            return
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        dx = (cur_xlim[1] - cur_xlim[0]) * 0.1
        dy = (cur_ylim[1] - cur_ylim[0]) * 0.1

        if event.key in ['w', 'W']:
            ax.set_ylim(cur_ylim[0] + dy, cur_ylim[1] + dy)
        elif event.key in ['s', 'S']:
            ax.set_ylim(cur_ylim[0] - dy, cur_ylim[1] - dy)
        elif event.key in ['a', 'A']:
            ax.set_xlim(cur_xlim[0] - dx, cur_xlim[1] - dx)
        elif event.key in ['d', 'D']:
            ax.set_xlim(cur_xlim[0] + dx, cur_xlim[1] + dx)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('key_press_event', on_key)

    def on_pick(evt):
        if hasattr(evt, 'mouseevent') and evt.mouseevent is not None:
            if evt.mouseevent.button != 1:  
                return
        
        if evt.artist == core_node:
            fig.canvas.mpl_disconnect(pick_cid) 
            show_sector_map(is_core=True)
            return 
        elif evt.artist == nodes_scatter:
            fig.canvas.mpl_disconnect(pick_cid) 
            show_sector_map(is_core=False)
            return

    pick_cid = fig.canvas.mpl_connect('pick_event', on_pick)
    fig.canvas.draw_idle()

def show_sector_map(is_core=False):
    """Draws the detailed local sector map with scroll zooming enabled."""
    fig.clf() 
    
    if is_core:
        num_stars = np.random.randint(600, 1200) 
    else:
        num_stars = np.random.randint(75000, 85000) 
        
    x = np.random.rand(num_stars) * 100
    y = np.random.rand(num_stars) * 100

    if is_core:
        distances = np.sqrt((x - 50)**2 + (y - 50)**2)
        valid_stars = distances > 25 
        x = x[valid_stars]
        y = y[valid_stars]
        num_stars = len(x)

    class_colors = np.array(['#5b82ff', '#9bb0ff', '#ffffff', '#fff4ea', '#ffe56f', '#ffb05b', '#ff4d4d'])
    class_weights = [0.001, 0.005, 0.015, 0.03, 0.08, 0.12, 0.749]
    star_classes = np.random.choice(len(class_colors), size=num_stars, p=class_weights)

    is_habitable = (np.random.rand(num_stars) < 0.006) & ((star_classes == 4) | (star_classes == 5))
    is_resource_rich = np.random.rand(num_stars) < 0.030
    is_exotic = np.random.rand(num_stars) < 0.001

    base_sizes = np.random.rand(num_stars) * 1.5 + 0.1
    base_colors = class_colors[star_classes]

    ax = fig.add_axes([0.22, 0.05, 0.75, 0.85])
    ax.set_facecolor('#050510')
    ax.axis('off')
    
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    title_text = f'Desolate Galactic Core ({num_stars:,} stars)' if is_core else f'Local Sector Map ({num_stars:,} stars)'
    title = ax.set_title(title_text + ' (Scroll to Zoom)', color='white', fontsize=16, fontfamily='monospace')

    scatter = ax.scatter(x, y, s=base_sizes, c=base_colors, alpha=0.6, edgecolors='none')

    if is_core:
        ax.scatter(50, 50, s=40000, color='#ff7700', alpha=0.08, zorder=8) 
        ax.scatter(50, 50, s=15000, color='#ffaa00', alpha=0.15, zorder=9) 
        ax.scatter(50, 50, s=5000, color='black', edgecolors='#ffffff', linewidths=2.0, zorder=10) 

    add_scroll_zoom(ax, max_range=300.0)

    ax_radio = plt.axes([0.02, 0.40, 0.17, 0.25], facecolor='#0a0a1a')
    radio = RadioButtons(ax_radio, ('Star Class', 'Habitable Systems', 'Resource Rich', 'Exotic Stars'), activecolor='#00ffcc')

    for label in radio.labels:
        label.set_color('white')
        label.set_fontsize(10)
        label.set_fontfamily('monospace')

    def set_mode(label):
        colors = base_colors.copy()
        sizes = base_sizes.copy()
        alphas = np.full(num_stars, 0.6)

        if label == 'Star Class':
            title.set_text(title_text + ' (Scroll to Zoom)')
        elif label == 'Habitable Systems':
            colors[:] = '#121222'; alphas[:] = 0.1; colors[is_habitable] = '#00ffcc'; sizes[is_habitable] = 14.0; alphas[is_habitable] = 1.0
        elif label == 'Resource Rich':
            colors[:] = '#121222'; alphas[:] = 0.1; colors[is_resource_rich] = '#ffd700'; sizes[is_resource_rich] = 8.0; alphas[is_resource_rich] = 0.9
        elif label == 'Exotic Stars':
            colors[:] = '#121222'; alphas[:] = 0.1; colors[is_exotic] = '#ff00ff'; sizes[is_exotic] = 20.0; alphas[is_exotic] = 1.0

        scatter.set_facecolors(colors)
        scatter.set_sizes(sizes)
        scatter.set_alpha(alphas)
        fig.canvas.draw_idle()

    radio.on_clicked(set_mode)

    ax_back = plt.axes([0.02, 0.85, 0.17, 0.05])
    btn_back = Button(ax_back, '< Return', color='#2a2a3a', hovercolor='#3a3a4a')
    btn_back.label.set_color('white')
    btn_back.label.set_fontfamily('monospace')
    btn_back.label.set_weight('bold')
    
    btn_back.on_clicked(show_galaxy_map)

    ui_elements['radio'] = radio
    ui_elements['btn_back'] = btn_back

    fig.canvas.draw_idle()

show_galaxy_map()
plt.show()