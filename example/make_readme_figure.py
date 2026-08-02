""" Generate the illustration figure of the README (example/readme_figure.svg):
	cartesian input -> polar (phi, r) map -> polar re-projection.
	Run from the repository root; requires numpy and matplotlib. """

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from CartToPolarDetector import convert

example_dir = os.path.dirname(os.path.abspath(__file__))
data_cart = np.load(os.path.join(example_dir, "input", "cartesian_detector.npz"))["example_cart_detector"]
N_x, N_y = data_cart.shape

N_phi, N_r = 120, 60
R_px_max = 0.5 * np.sqrt(N_x**2 + N_y**2)
data_polar = convert("readme_figure", data_cart, N_phi=N_phi, N_r=N_r)

# flux -> surface brightness (divide by polar pixel area) for a fair visual comparison
dphi = 2. * np.pi / N_phi
r_edges = np.linspace(0., R_px_max, N_r + 1)
areas = 0.5 * dphi * (r_edges[1:]**2 - r_edges[:-1]**2)
data_polar_sb = data_polar / areas[None, :]

phi_edges = np.linspace(0., 2. * np.pi, N_phi + 1)

fig = plt.figure(figsize=(12, 4.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.15, 1], wspace=0.35)
vmax = np.max(data_cart)

# 1. cartesian input
ax1 = fig.add_subplot(gs[0])
ax1.imshow(data_cart.T, origin="lower", cmap="inferno", vmin=0, vmax=vmax,
	extent=[-N_x/2, N_x/2, -N_y/2, N_y/2])
ax1.set_title("cartesian detector (input)")
ax1.set_xlabel("$x$ [px]")
ax1.set_ylabel("$y$ [px]")

# 2. polar detector as (phi, r) map
ax2 = fig.add_subplot(gs[1])
im2 = ax2.pcolormesh(np.degrees(phi_edges), r_edges, data_polar_sb.T,
	cmap="inferno", vmin=0, vmax=vmax, rasterized=True)
ax2.set_title(r"polar detector ($\varphi$, $r$) map")
ax2.set_xlabel(r"$\varphi$ [deg]")
ax2.set_ylabel("$r$ [px]")
ax2.set_xticks([0, 90, 180, 270, 360])
cb = fig.colorbar(im2, ax=ax2, pad=0.02)
cb.set_label("surface brightness")

# 3. polar re-projection
ax3 = fig.add_subplot(gs[2], projection="polar")
ax3.set_theta_zero_location("N")	# phi = 0 points in +y direction
ax3.set_theta_direction(-1)			# phi increases towards +x
ax3.pcolormesh(phi_edges, r_edges, data_polar_sb.T,
	cmap="inferno", vmin=0, vmax=vmax, rasterized=True)
ax3.set_title("polar detector (re-projected)")
ax3.set_yticklabels([])
ax3.grid(alpha=0.25)

fig.savefig(os.path.join(example_dir, "readme_figure.svg"), dpi=72, bbox_inches="tight")
print("flux conserved:", np.isclose(data_cart.sum(), data_polar.sum()))
