""" Result-value tests for CartToPolarDetector.

Opposed to the invariant tests, these tests check the actual numbers that
convert() produces, by comparing them with independently obtained references:
- exact analytical intersection areas for simple geometries
  (quadrants, circles inside a pixel, a circular segment),
- a brute-force supersampling integration of the very same conversion,
  which shares no code with the tested implementation.
"""

import os

import numpy as np
import pytest

from CartToPolarDetector import convert

HERE = os.path.dirname(__file__)
EXAMPLE_INPUT = os.path.join(HERE, "..", "example", "input", "cartesian_detector.npz")

pix2 = 2. * np.pi


def polar_reference(data_cart, N_phi, N_r, center_shift, R_px_max, s=100):
	""" independent brute-force reference: every cartesian pixel is subdivided
		into s*s subpixels of equal flux, and each subpixel is assigned to the
		polar cell that contains its center """
	N_x, N_y = data_cart.shape
	center = np.array([N_x / 2., N_y / 2.]) + np.array(center_shift)
	dphi, dr = pix2 / N_phi, R_px_max / N_r

	x = (np.arange(N_x * s) + 0.5) / s - center[0]
	y = (np.arange(N_y * s) + 0.5) / s - center[1]
	X, Y = np.meshgrid(x, y, indexing="ij")
	phi = np.arctan2(X, Y) % pix2
	r = np.hypot(X, Y)
	weights = np.repeat(np.repeat(data_cart, s, axis=0), s, axis=1) / (s * s)

	mask = r < R_px_max
	n_phi = np.minimum((phi[mask] / dphi).astype(int), N_phi - 1)
	n_r = np.minimum((r[mask] / dr).astype(int), N_r - 1)
	reference = np.zeros((N_phi, N_r))
	np.add.at(reference, (n_phi, n_r), weights[mask])
	return reference


def test_quadrant_mapping_exact():
	""" for a 2x2 detector with N_phi=4 and one radial bin, every cartesian
		pixel lies exactly in one polar quadrant (phi=0 points in +y direction,
		phi=pi/2 in +x direction) """
	data_cart = np.array([[1., 2.], [3., 4.]])
	data_polar = convert("test", data_cart, N_phi=4, N_r=1)
	assert data_polar.shape == (4, 1)
	# quadrant 0 (+x,+y): pixel [1,1]; quadrant 1 (+x,-y): pixel [1,0];
	# quadrant 2 (-x,-y): pixel [0,0]; quadrant 3 (-x,+y): pixel [0,1]
	assert np.allclose(data_polar[:, 0], [4., 3., 1., 2.], rtol=0., atol=1e-12)


def test_circle_inside_central_pixel_exact():
	""" a single pixel with the polar center in its middle: the innermost ring
		is a full circle inside the pixel, so each of its cells has the exact
		area pi*dr^2/4; the outer ring holds the rest of the quadrant """
	value = 3.7
	data_cart = np.array([[value]])
	N_phi, N_r = 4, 2
	R_px_max = 0.5 * np.sqrt(2.)			# default for a 1x1 detector
	dr = R_px_max / N_r						# dr < 0.5, circle fully inside
	data_polar = convert("test", data_cart, N_phi=N_phi, N_r=N_r)

	assert np.allclose(data_polar[:, 0], value * np.pi * dr * dr / 4., rtol=1e-12)
	assert np.allclose(data_polar[:, 1], value * (0.25 - np.pi * dr * dr / 4.), rtol=1e-12)


def test_circular_segment_exact():
	""" a single bright pixel on the +y axis cut by a circle that crosses only
		its bottom edge: the inner part is exactly the circular segment
		A = r0^2*arccos(h/r0) - h*sqrt(r0^2 - h^2) above the chord at height h """
	value = 3.7
	data_cart = np.zeros((1, 5))
	data_cart[0, 4] = value					# pixel center at c = (0, 2)
	N_phi, N_r = 4, 2
	r0, h = 1.55, 1.5						# circle radius and bottom edge height
	R_px_max = N_r * r0						# makes r0 the border of ring 0
	data_polar = convert("test", data_cart, N_phi=N_phi, N_r=N_r, R_px_max_in=R_px_max)

	segment = r0 * r0 * np.arccos(h / r0) - h * np.sqrt(r0 * r0 - h * h)
	assert np.allclose(data_polar[:, 0].sum(), value * segment, rtol=1e-12)
	assert np.allclose(data_polar[:, 1].sum(), value * (1. - segment), rtol=1e-12)
	# the pixel is symmetric about the +y axis (phi=0), so the flux is split
	# equally between the first and the last phi bin
	assert np.allclose(data_polar[0, :], data_polar[-1, :], rtol=1e-12)


@pytest.mark.parametrize("N_phi,N_r,center_shift,R_px_max", [
	(8, 6, [0., 0.], 15.),
	(12, 5, [1.5, -3.], 10.),
])
def test_against_supersampling_reference(N_phi, N_r, center_shift, R_px_max):
	""" the full result on the example data agrees with an independent
		brute-force integration of the same polar grid """
	data_cart = np.load(EXAMPLE_INPUT)["example_cart_detector"]
	data_polar = convert("test", data_cart, N_phi=N_phi, N_r=N_r,
		center_shift=center_shift, R_px_max_in=R_px_max)
	reference = polar_reference(data_cart, N_phi, N_r, center_shift, R_px_max, s=100)

	# the reference resolves polar cell borders only up to the subpixel size,
	# hence the moderate tolerance
	assert np.max(np.abs(data_polar - reference)) < 5e-3 * np.max(np.abs(data_polar))
