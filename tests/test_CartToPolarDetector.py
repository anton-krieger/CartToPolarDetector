""" Tests for CartToPolarDetector.

They lock in the invariants of the conversion method (see Appendix B of
Krieger & Wolf 2022, A&A, 662, A99):
- the overall flux is conserved,
- for a uniform input, every polar pixel value equals its geometrical area,
- rotating the input by 90 degrees shifts the output by N_phi/4,
- the stored example result is reproduced exactly,
- input files are read back correctly and invalid arguments are rejected.

Run from the repository root with: pytest
"""

import os

import numpy as np
import pytest

from CartToPolarDetector import convert, read_input

HERE = os.path.dirname(__file__)
EXAMPLE_INPUT = os.path.join(HERE, "..", "example", "input", "cartesian_detector.npz")
EXAMPLE_RESULT = os.path.join(HERE, "..", "example", "results", "wave_example3")


@pytest.fixture
def example_data():
	""" the cartesian example detector shipped in example/input """
	return np.load(EXAMPLE_INPUT)["example_cart_detector"]


def test_flux_conservation_example_data(example_data):
	data_polar = convert("test", example_data, N_phi=96, N_r=96)
	assert abs(data_polar.sum() - example_data.sum()) < 1e-9 * np.abs(example_data).sum()


def test_flux_conservation_uniform_shifted_center():
	data_cart = np.ones((17, 23))
	data_polar = convert("test", data_cart, N_phi=16, N_r=13, center_shift=[0.3, -0.7])
	assert abs(data_polar.sum() - data_cart.sum()) < 1e-9 * data_cart.sum()
	assert data_polar.min() > -1e-12


def test_uniform_input_gives_polar_pixel_areas():
	""" for a uniform input of ones, the value of every polar pixel that lies
		fully inside the detector equals its geometrical area 0.5*dphi*(r_out^2 - r_in^2) """
	N = 20
	N_phi, N_r = 16, 10
	data_polar = convert("test", np.ones((N, N)), N_phi=N_phi, N_r=N_r)

	R_px_max = 0.5 * np.sqrt(2.) * N
	dphi = 2. * np.pi / N_phi
	dr = R_px_max / N_r
	for n_r in range(N_r):
		r_in, r_out = n_r * dr, (n_r + 1) * dr
		if r_out <= N / 2.:  # ring fully inside the detector
			expected = 0.5 * dphi * (r_out * r_out - r_in * r_in)
			assert np.allclose(data_polar[:, n_r], expected, rtol=1e-9)


def test_rotation_symmetry(example_data):
	""" rotating the input by 90 degrees shifts the output by N_phi/4 """
	N_phi, N_r = 32, 20
	p1 = convert("test", example_data, N_phi=N_phi, N_r=N_r)
	p2 = convert("test", np.rot90(example_data, k=-1), N_phi=N_phi, N_r=N_r)
	assert np.allclose(np.roll(p1, N_phi // 4, axis=0), p2)


def test_reproduces_shipped_example_result(example_data):
	""" regression test: the parameters stored in example/results/wave_example3/input.dat
		reproduce the stored output.npz """
	inp = read_input(os.path.join(EXAMPLE_RESULT, "input.dat"))
	reference = np.load(os.path.join(EXAMPLE_RESULT, "output.npz"))["polar"]
	data_polar = convert("test", example_data, N_phi=inp["N_phi"], N_r=inp["N_r"],
		center_shift=inp["center_shift"], R_px_max_in=inp["R_px_max_in"])
	assert data_polar.shape == reference.shape
	assert np.allclose(data_polar, reference)


def test_save_and_read_input_roundtrip(example_data, tmp_path):
	convert("test", example_data, N_phi=8, N_r=5, center_shift=[0.25, -0.5],
		results_folder_path=str(tmp_path), save_name="roundtrip", progress_bar=False)
	saved = tmp_path / "roundtrip"

	inp = read_input(saved / "input.dat")
	assert inp["N_phi"] == 8
	assert inp["N_r"] == 5
	assert inp["center_shift"] == [0.25, -0.5]
	assert isinstance(inp["dr"], float)
	assert isinstance(inp["dphi"], float)
	assert inp["R_px_max_in"] is None
	assert inp["progress_bar"] is False

	assert len((saved / "phi_border_values.dat").read_text().split()) == 8 + 1
	assert len((saved / "r_border_values.dat").read_text().split()) == 5 + 1


def test_invalid_arguments_raise(example_data, tmp_path):
	with pytest.raises(ValueError):
		convert("test", example_data, N_phi=30, N_r=10)  # N_phi not a multiple of 4
	with pytest.raises(ValueError):
		convert("test", example_data, N_phi=8, N_r=0)  # N_r not positive
	with pytest.raises(ValueError):
		convert("test", example_data, N_phi=8, N_r=5, results_folder_path=str(tmp_path / "missing"))
