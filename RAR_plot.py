import numpy as np
import matplotlib.pyplot as plt

plot_acc = False

x = np.array([1, 2, 3, 4, 5], dtype=float)
f = np.array([3, 5, 3.5, 2.5, 1], dtype=float)
g = np.array([5, 4, 3, 2, 1], dtype=float)

selected_idx = None
selected_curve = None  # 'f' or 'g'

# --- Ranking with tie handling (average ranks) ---
def rankdata(a):
    temp = np.argsort(a)
    ranks = np.empty_like(temp, dtype=float)
    ranks[temp] = np.arange(len(a))

    # Handle ties (average ranks)
    unique, inverse, counts = np.unique(a, return_inverse=True, return_counts=True)
    for u, c in zip(unique, counts):
        if c > 1:
            idx = np.where(a == u)[0]
            ranks[idx] = ranks[idx].mean()
    return ranks + 1  # 1-based ranks

def spearman_rho(f, g):
    rf = rankdata(f)
    rg = rankdata(g)
    return np.corrcoef(rf, rg)[0, 1]

def pearson_rho(f, g):
    return np.corrcoef(f, g)[0, 1]

def compute_sorted(f, g):
    idx = np.argsort(g)
    return g[idx], f[idx]

# --- Setup figure ---
fig, axs = plt.subplots(3, 1, figsize=(6, 8))
ax_f, ax_g, ax_fg = axs

# f(x)
line_f, = ax_f.plot(x, f, 'o-', picker=5)
ax_f.set_title("f(x) vs x")

# g(x)
line_g, = ax_g.plot(x, g, 'o-', picker=5)
ax_g.set_title("g(x) vs x")

# f vs g (sorted)
if plot_acc:
    g_sorted, f_sorted = compute_sorted(f**2/x, g**2/x)
    line_fg, = ax_fg.plot(g_sorted, f_sorted, 'o-')
    ax_fg.set_title("f^2/x vs g^2/x (sorted by g^2/x)")
    ax_fg.set_xlabel("g^2/x")
    ax_fg.set_ylabel("f^2/x")
else:
    g_sorted, f_sorted = compute_sorted(f, g)
    line_fg, = ax_fg.plot(g_sorted, f_sorted, 'o-')
    ax_fg.set_title("f vs g (sorted by g)")
    ax_fg.set_xlabel("g")
    ax_fg.set_ylabel("f")

# Spearman text
rho_text = ax_fg.text(
    0.05, 0.95, "", transform=ax_fg.transAxes,
    verticalalignment='top', fontsize=12
)

def update_plots():
    line_f.set_ydata(f)
    line_g.set_ydata(g)

    if plot_acc: g_sorted, f_sorted = compute_sorted(f**2/x, g**2/x)
    else: g_sorted, f_sorted = compute_sorted(f, g)
    line_fg.set_xdata(g_sorted)
    line_fg.set_ydata(f_sorted)

    # Compute Spearman and Pearson rho
    rho_spearman = spearman_rho(f, g)
    rho_pearson = pearson_rho(f, g)

    # Update text + color cue
    if np.isclose(rho_spearman, 1):
        color = "green"
    elif np.isclose(rho_pearson, -1):
        color = "blue"
    else:
        color = "red"

    rho_text.set_text(f"Spearman ρ = {rho_spearman:.3f}\nPearson ρ = {rho_pearson:.3f}")
    rho_text.set_color(color)

    for ax in axs:
        ax.relim()
        ax.autoscale_view()

    fig.canvas.draw_idle()

def on_pick(event):
    global selected_idx, selected_curve

    artist = event.artist
    mouse_event = event.mouseevent

    if artist == line_f:
        selected_curve = 'f'
    elif artist == line_g:
        selected_curve = 'g'
    else:
        return

    xdata = artist.get_xdata()
    ydata = artist.get_ydata()

    distances = np.hypot(xdata - mouse_event.xdata, ydata - mouse_event.ydata)
    selected_idx = np.argmin(distances)

def on_motion(event):
    global selected_idx, selected_curve

    if selected_idx is None or event.inaxes is None:
        return

    if selected_curve == 'f' and event.inaxes == ax_f:
        f[selected_idx] = event.ydata
    elif selected_curve == 'g' and event.inaxes == ax_g:
        g[selected_idx] = event.ydata
    else:
        return

    update_plots()

def on_release(event):
    global selected_idx, selected_curve
    selected_idx = None
    selected_curve = None

fig.canvas.mpl_connect('pick_event', on_pick)
fig.canvas.mpl_connect('motion_notify_event', on_motion)
fig.canvas.mpl_connect('button_release_event', on_release)

update_plots()
plt.tight_layout()
plt.show()
