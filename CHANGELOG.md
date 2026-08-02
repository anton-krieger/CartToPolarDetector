# Changelog

## 0.1.0 (unreleased)

### Fixed
- `NameError` in error branches of `get_ray_pixel_intersect_M` (b/c/d cases referenced the wrong variable).
- Missing parentheses on `sys.exit` calls, which made them no-ops.
- Removed a duplicate dead `elif` branch in `get_area_summed_TR`.
- `read_input` now correctly parses `progress_bar` as a boolean and `R_px_max_in` as a float or `None` (previously both were returned as strings).
- Errors now raise proper exceptions (`ValueError` for invalid user input, `RuntimeError` for internal invariants) instead of printing and exiting.

### Changed
- Renamed misspelled functions: `get_rotated_pixel_contributionn` -> `get_rotated_pixel_contribution` and `summed_to_indivudual_contribution` -> `summed_to_individual_contribution`. The old names remain available as aliases.
- Deduplicated the four per-edge code blocks in `get_ring_pos`, `get_circle_pixel_intersect`, and `get_area_summed_M` into loops. Results are bit-identical.
- Performance: replaced `np.pad`-based finite differences with slice arithmetic and scalar `numpy` math calls with the `math` module (about 2x faster overall). Results are bit-identical.
- Packaging modernized: metadata moved to `pyproject.toml` (setuptools backend), `setup.py` reduced to a compatibility shim, version bumped to 0.1.0, `requires-python >= 3.9`.

### Added
- Test suite (`tests/`) with invariant tests (flux conservation, rotation symmetry, round trips, input validation) and result-value tests (exact analytical areas and an independent brute-force supersampling reference).
- GitHub Actions workflow running the tests on Python 3.9, 3.11, and 3.13.
- Rewritten `README.md` with installation, quick start, examples, tests, and citation sections.
- README illustration figure (`example/readme_figure.svg`) and its generation script (`example/make_readme_figure.py`).
- Module docstring referencing the method description (Krieger & Wolf 2022, A&A, 662, A99, Appendix B).

## 0.0.1 (2022-01-17)

- Initial release.
