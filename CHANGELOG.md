# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-07-23

### Added
- Full Python port of the Pan-Tompkins QRS detection algorithm
- `pan_tompkins(ecg, fs, gr)` as the main detection function
- Support for arbitrary sampling rates via adaptive filter kernel interpolation
- Separate LP + HP filtering path for 200 Hz signals; bandpass for all other rates
- Zero-phase filtering using `scipy.signal.filtfilt` throughout
- Adaptive dual-threshold decision logic with T-wave rejection
- Search-back procedure for missed QRS complexes
- Optional Matplotlib visualisation of all intermediate processing stages
- Test suite: 16 tests across 5 classes (`tests/test_detector.py`)
- `pyproject.toml` with Hatch build backend for PyPI compatibility
- GitHub Actions CI: automated testing on Python 3.9–3.12 across Ubuntu, macOS, Windows
- GitHub Actions workflow for auto-compiling JOSS paper PDF on push
- Quickstart tutorial (`docs/quickstart.md`)
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `CHANGELOG.md`, `CONTRIBUTORS`
- BSD 3-Clause license with full attribution to original MATLAB author Hooman Sedghamiz
