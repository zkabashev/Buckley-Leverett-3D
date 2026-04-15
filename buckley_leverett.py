import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from matplotlib.collections import LineCollection
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  ROCK / FLUID MODEL
# ─────────────────────────────────────────────────────────────────────────────
S_WC  = 0.20          # connate water saturation
S_OR  = 0.20          # residual oil saturation
MU_O  = 1.0           # oil viscosity  [cP]  – fixed
N_W   = 2.0           # Corey exponent – water
N_O   = 2.0           # Corey exponent – oil
U_TOT = 1.0           # total Darcy velocity (normalised)

def normalised_sat(sw):
    """Effective (normalised) water saturation."""
    return np.clip((sw - S_WC) / (1.0 - S_WC - S_OR), 0.0, 1.0)

def rel_perms(sw):
    sn = normalised_sat(sw)
    krw = sn ** N_W
    kro = (1.0 - sn) ** N_O
    return krw, kro

def fractional_flow(sw, mu_w):
    """Water fractional flow  f_w(S_w ; μ_w)."""
    krw, kro = rel_perms(sw)
    lam_w = krw / mu_w
    lam_o = kro / MU_O
    denom = lam_w + lam_o
    # avoid division by zero at endpoints
    fw = np.where(denom > 1e-12, lam_w / denom, 0.0)
    return fw

# ─────────────────────────────────────────────────────────────────────────────
# 2.  NUMERICAL SOLVER  (explicit upwind finite-difference)
# ─────────────────────────────────────────────────────────────────────────────
NX   = 200          # spatial cells
NT   = 300          # time steps stored
XMAX = 1.0          # dimensionless reservoir length
TMAX = 1.0          # dimensionless pore-volumes injected

def run_simulation(mu_w: float):
    """Return saturation array S[nx, nt] for given μ_w."""
    dx = XMAX / NX
    # CFL-limited time step  (max df/dS ≤ ~ 4 for typical Corey)
    max_dfw = 4.0
    dt_cfl  = 0.45 * dx / (U_TOT * max_dfw)

    # choose stored snapshots evenly spaced in time
    t_out   = np.linspace(0.0, TMAX, NT)
    dt_store = t_out[1] - t_out[0]
    n_sub   = max(1, int(np.ceil(dt_store / dt_cfl)))
    dt      = dt_store / n_sub          # actual integration dt

    x = np.linspace(0.5 * dx, XMAX - 0.5 * dx, NX)
    S = np.full(NX, S_WC)              # initial condition: connate water

    S_out = np.zeros((NX, NT))
    S_out[:, 0] = S.copy()

    for k in range(1, NT):
        for _ in range(n_sub):
            fw  = fractional_flow(S, mu_w)
            # upwind flux at cell faces (injection from left at fw = 1)
            flux_left  = np.empty(NX)
            flux_left[0]  = fractional_flow(np.array([1.0 - S_OR]), mu_w)[0]
            flux_left[1:] = fw[:-1]          # upwind: flow left → right
            dS = -(U_TOT / dx) * (fw - flux_left) * dt
            S  = np.clip(S + dS, S_WC, 1.0 - S_OR)
        S_out[:, k] = S.copy()

    return x, t_out, S_out

# ─────────────────────────────────────────────────────────────────────────────
# 3.  FIGURE LAYOUT
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 9), facecolor="#111111")
fig.suptitle("1D Buckley–Leverett Waterflood  ·  Saturation Evolution",
             color="white", fontsize=14, fontweight="bold", y=0.98)

gs = gridspec.GridSpec(
    3, 3,
    left=0.05, right=0.97, top=0.93, bottom=0.22,
    wspace=0.35, hspace=0.35,
    height_ratios=[0.05, 1.0, 0.7]
)

# Equation label row (row 0 spans all columns)
ax_eq = fig.add_subplot(gs[0, :])
ax_eq.axis("off")
ax_eq.text(0.5, 0.5,
           r"$\dfrac{\partial S_w}{\partial t} + u\,\dfrac{\partial f_w}{\partial x} = 0$"
           r"$\quad f_w = \dfrac{k_{rw}/\mu_w}{k_{rw}/\mu_w + k_{ro}/\mu_o}$",
           ha="center", va="center", color="lightyellow",
           fontsize=13, transform=ax_eq.transAxes)

# 3D surface (row 1, col 0:2)
ax3d = fig.add_subplot(gs[1, :2], projection="3d")
ax3d.set_facecolor("#1a1a1a")
for a in [ax3d.xaxis, ax3d.yaxis, ax3d.zaxis]:
    a.pane.fill = False
    a.line.set_color("#444444")

# Time-slice 2-D profile (row 2, col 0:2)  ← under the 3D plot
ax_t = fig.add_subplot(gs[2, :2])
ax_t.set_facecolor("#1a1a1a")
ax_t.tick_params(colors="white"); ax_t.xaxis.label.set_color("white")
ax_t.yaxis.label.set_color("white")
for sp in ax_t.spines.values(): sp.set_color("#555555")

# Space-slice history (row 1, col 2)  ← beside the 3D plot
ax_x = fig.add_subplot(gs[1, 2])
ax_x.set_facecolor("#1a1a1a")
ax_x.tick_params(colors="white"); ax_x.xaxis.label.set_color("white")
ax_x.yaxis.label.set_color("white")
for sp in ax_x.spines.values(): sp.set_color("#555555")

# Fractional-flow curve (row 2, col 2)
ax_fw = fig.add_subplot(gs[2, 2])
ax_fw.set_facecolor("#1a1a1a")
ax_fw.tick_params(colors="white"); ax_fw.xaxis.label.set_color("white")
ax_fw.yaxis.label.set_color("white")
for sp in ax_fw.spines.values(): sp.set_color("#555555")

# ─────────────────────────────────────────────────────────────────────────────
# 4.  SLIDER AXES  (below the plots)
# ─────────────────────────────────────────────────────────────────────────────
slider_color  = "#2a2a2a"
handle_color  = "#00ccff"

def make_slider(rect, label, vmin, vmax, vinit, fmt="%0.3f"):
    ax_s = plt.axes(rect, facecolor=slider_color)
    sl   = Slider(ax_s, label, vmin, vmax, valinit=vinit, valfmt=fmt,
                  color=handle_color)
    sl.label.set_color("white"); sl.valtext.set_color("white")
    return sl

sl_t   = make_slider([0.08, 0.14, 0.35, 0.025],
                      "Time  (PVI)",   0.0,  1.0, 0.5)
sl_x   = make_slider([0.08, 0.10, 0.35, 0.025],
                      "Space (x/L)",   0.0,  1.0, 0.5)
sl_mu  = make_slider([0.08, 0.06, 0.35, 0.025],
                      "μ_w  (cP)",      0.2, 10.0, 1.0, fmt="%0.2f")

# Info text
info_ax = plt.axes([0.55, 0.04, 0.40, 0.12], facecolor="#1a1a1a")
info_ax.axis("off")
info_text = info_ax.text(
    0.05, 0.95,
    "Drag sliders to explore the saturation field.\n"
    "• Time slider  → vertical plane (S_w vs x at fixed t)\n"
    "• Space slider → lateral plane  (S_w vs t at fixed x)\n"
    "• μ_w slider   → reruns simulation with new viscosity",
    color="#aaaaaa", fontsize=8.5, va="top", wrap=True,
    transform=info_ax.transAxes
)

# ─────────────────────────────────────────────────────────────────────────────
# 5.  COLOURMAP helper
# ─────────────────────────────────────────────────────────────────────────────
CMAP = plt.cm.plasma        # similar hot-to-cool look

def sat_color(s_val):
    """Map saturation [S_WC, 1-S_OR] to RGBA."""
    lo, hi = S_WC, 1.0 - S_OR
    return CMAP((s_val - lo) / (hi - lo))

# ─────────────────────────────────────────────────────────────────────────────
# 6.  SHOCK CONSTRUCTION
# ─────────────────────────────────────────────────────────────────────────────
def _shock_saturation(mu_w, n_pts=500):
    """Welge tangent construction – returns shock-front saturation."""
    sw_arr = np.linspace(S_WC + 1e-4, 1.0 - S_OR - 1e-4, n_pts)
    fw_arr = fractional_flow(sw_arr, mu_w)
    # slope from (S_wc, 0)
    slope  = (fw_arr - 0.0) / (sw_arr - S_WC)
    return sw_arr[np.argmax(slope)]


# ─────────────────────────────────────────────────────────────────────────────
# 7.  DRAW / REDRAW FUNCTION
# ─────────────────────────────────────────────────────────────────────────────
_cache = {}   # store last simulation results

def draw(mu_w, t_val, x_val):
    global _cache

    # re-run only when μ_w changes
    if _cache.get("mu_w") != mu_w:
        x, t, S = run_simulation(mu_w)
        _cache = dict(mu_w=mu_w, x=x, t=t, S=S)
    else:
        x, t, S = _cache["x"], _cache["t"], _cache["S"]

    t_idx = int(np.argmin(np.abs(t - t_val)))
    x_idx = int(np.argmin(np.abs(x - x_val)))

    T_grid, X_grid = np.meshgrid(t, x)      # shapes (NX, NT)

    # ── 3-D surface ──────────────────────────────────────────────────────────
    ax3d.cla()
    ax3d.set_facecolor("#1a1a1a")
    for a in [ax3d.xaxis, ax3d.yaxis, ax3d.zaxis]:
        a.pane.fill = False
        a.line.set_color("#444444")

    # Subsample for speed
    step_x, step_t = max(1, NX // 60), max(1, NT // 60)
    Xs = X_grid[::step_x, ::step_t]
    Ts = T_grid[::step_x, ::step_t]
    Ss = S[::step_x, ::step_t]

    fcolors = CMAP((Ss - S_WC) / (1 - S_WC - S_OR))
    surf = ax3d.plot_surface(Xs, Ts, Ss, facecolors=fcolors,
                             linewidth=0, antialiased=True, alpha=0.82)

    # Wireframe overlay (thin)
    ax3d.plot_wireframe(Xs, Ts, Ss, color="white", linewidth=0.15, alpha=0.25)

    # Time-slice plane
    s_tslice = S[:, t_idx]
    ax3d.plot(x, np.full_like(x, t[t_idx]), s_tslice,
              color="#00ffff", linewidth=2.0, zorder=10)
    ax3d.plot_surface(
        np.column_stack([x, x]).reshape(NX, 2),
        np.full((NX, 2), t[t_idx]),
        np.column_stack([np.zeros(NX), s_tslice]).reshape(NX, 2),
        color="#00ffff", alpha=0.10
    )

    # Space-slice plane
    s_xslice = S[x_idx, :]
    ax3d.plot(np.full_like(t, x[x_idx]), t, s_xslice,
              color="#ffaa00", linewidth=2.0, zorder=10)
    ax3d.plot_surface(
        np.full((2, NT), x[x_idx]),
        np.row_stack([t, t]),
        np.row_stack([np.zeros(NT), s_xslice]),
        color="#ffaa00", alpha=0.10
    )

    ax3d.set_xlabel("x / L", color="white", labelpad=6, fontsize=9)
    ax3d.set_ylabel("Time (PVI)", color="white", labelpad=6, fontsize=9)
    ax3d.set_zlabel("Water Saturation  Sw", color="white", labelpad=6, fontsize=9)
    ax3d.set_title("Saturation Surface  S_w(x, t)", color="white",
                   fontsize=10, pad=4)
    ax3d.tick_params(colors="white", labelsize=7)
    ax3d.set_zlim(0, 1)
    ax3d.view_init(elev=28, azim=-55)

    # ── Time-slice panel ─────────────────────────────────────────────────────
    ax_t.cla(); ax_t.set_facecolor("#1a1a1a")
    # gradient colouring along saturation
    points = np.array([x, s_tslice]).T.reshape(-1, 1, 2)
    segs   = np.concatenate([points[:-1], points[1:]], axis=1)
    norm   = plt.Normalize(S_WC, 1 - S_OR)
    lc     = LineCollection(segs, cmap=CMAP, norm=norm, linewidth=2.5)
    lc.set_array(s_tslice[:-1])
    ax_t.add_collection(lc)
    ax_t.fill_between(x, S_WC, s_tslice, alpha=0.25, color="#00ffff")
    ax_t.axhline(S_WC, color="#555555", linewidth=0.8, linestyle="--")
    ax_t.set_xlim(0, XMAX); ax_t.set_ylim(0, 1)
    ax_t.set_xlabel("x / L", fontsize=9); ax_t.set_ylabel("Sw", fontsize=9)
    ax_t.set_title(f"Profile at  t = {t[t_idx]:.3f} PVI",
                   color="#00ffff", fontsize=9)
    ax_t.tick_params(colors="white", labelsize=7)
    for sp in ax_t.spines.values(): sp.set_color("#555555")

    # Shock front indicator – Rankine-Hugoniot: v_shock = Δf_w / ΔS_w
    sw_shock   = _shock_saturation(mu_w)
    fw_shock   = fractional_flow(np.array([sw_shock]), mu_w)[0]
    v_shock    = (fw_shock - 0.0) / (sw_shock - S_WC)   # speed (dim-less)
    x_front    = v_shock * t[t_idx]                      # position at this t
    x_front    = min(x_front, XMAX)                      # clamp to domain
    ax_t.axvline(x_front, color="red", linewidth=1.4,
                 linestyle=":", label=f"Front  x={x_front:.3f}")
    ax_t.legend(fontsize=7, facecolor="#222222", labelcolor="white",
                edgecolor="#555555")

    # ── Space-slice panel ────────────────────────────────────────────────────
    ax_x.cla(); ax_x.set_facecolor("#1a1a1a")
    points2 = np.array([t, s_xslice]).T.reshape(-1, 1, 2)
    segs2   = np.concatenate([points2[:-1], points2[1:]], axis=1)
    lc2     = LineCollection(segs2, cmap=CMAP, norm=norm, linewidth=2.5)
    lc2.set_array(s_xslice[:-1])
    ax_x.add_collection(lc2)
    ax_x.fill_between(t, S_WC, s_xslice, alpha=0.25, color="#ffaa00")
    ax_x.axvline(t[t_idx], color="#00ffff", linewidth=1.0,
                 linestyle="--", label=f"t = {t[t_idx]:.3f}")
    ax_x.set_xlim(0, TMAX); ax_x.set_ylim(0, 1)
    ax_x.set_xlabel("Time (PVI)", fontsize=9)
    ax_x.set_ylabel("Sw", fontsize=9)
    ax_x.set_title(f"History at  x = {x[x_idx]:.3f} L",
                   color="#ffaa00", fontsize=9)
    ax_x.tick_params(colors="white", labelsize=7)
    ax_x.legend(fontsize=7, facecolor="#222222", labelcolor="white",
                edgecolor="#555555")
    for sp in ax_x.spines.values(): sp.set_color("#555555")

    # ── Fractional-flow panel ────────────────────────────────────────────────
    ax_fw.cla(); ax_fw.set_facecolor("#1a1a1a")
    sw_arr = np.linspace(S_WC, 1 - S_OR, 300)
    fw_arr = fractional_flow(sw_arr, mu_w)
    ax_fw.plot(sw_arr, fw_arr, color="#99ff66", linewidth=2)

    # Welge tangent
    sw_sh = _shock_saturation(mu_w)
    fw_sh = fractional_flow(np.array([sw_sh]), mu_w)[0]
    ax_fw.plot([S_WC, sw_sh], [0, fw_sh], "--", color="#ff4444",
               linewidth=1.5, label=f"Shock  Sw={sw_sh:.2f}")
    ax_fw.scatter([sw_sh], [fw_sh], color="red", zorder=5, s=30)
    ax_fw.set_xlabel("Sw", fontsize=9); ax_fw.set_ylabel("f_w", fontsize=9)
    ax_fw.set_title(f"Fractional Flow  (μ_w={mu_w:.2f} cP)",
                    color="#99ff66", fontsize=9)
    ax_fw.set_xlim(0, 1); ax_fw.set_ylim(0, 1)
    ax_fw.tick_params(colors="white", labelsize=7)
    ax_fw.legend(fontsize=7, facecolor="#222222", labelcolor="white",
                 edgecolor="#555555")
    for sp in ax_fw.spines.values(): sp.set_color("#555555")

    fig.canvas.draw_idle()


# ─────────────────────────────────────────────────────────────────────────────
# 8.  SLIDER CALLBACKS
# ─────────────────────────────────────────────────────────────────────────────
def on_change(_):
    draw(sl_mu.val, sl_t.val, sl_x.val)

sl_t.on_changed(on_change)
sl_x.on_changed(on_change)
sl_mu.on_changed(on_change)

# ─────────────────────────────────────────────────────────────────────────────
# 9.  INITIAL DRAW
# ─────────────────────────────────────────────────────────────────────────────
draw(mu_w=1.0, t_val=0.5, x_val=0.5)

plt.show()
