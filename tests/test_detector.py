"""
Tests for pantompkins.pan_tompkins
"""

import numpy as np
import pytest
from pantompkins import pan_tompkins


def synthetic_ecg(fs: float, duration: float = 10.0, bpm: float = 70.0, noise: float = 0.05):
    """Generate a synthetic ECG-like signal with narrow Gaussian QRS complexes."""
    t = np.linspace(0, duration, int(duration * fs))
    ecg = np.zeros_like(t)
    beat_interval = 60.0 / bpm
    beat_times = np.arange(beat_interval / 2, duration, beat_interval)
    for bt in beat_times:
        ecg += np.exp(-((t - bt) ** 2) / (2 * (0.02) ** 2))
    rng = np.random.default_rng(42)
    ecg += noise * rng.standard_normal(len(t))
    return ecg, len(beat_times)


# ── Basic detection ──────────────────────────────────────────────────────────

class TestBasicDetection:

    def test_returns_three_values(self):
        ecg, _ = synthetic_ecg(200)
        result = pan_tompkins(ecg, 200, gr=False)
        assert len(result) == 3

    def test_output_shapes_consistent(self):
        ecg, _ = synthetic_ecg(200)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert amp.shape == idx.shape

    def test_detects_reasonable_beat_count_200hz(self):
        ecg, expected = synthetic_ecg(200)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert abs(len(idx) - expected) <= 2, f"Expected ~{expected} beats, got {len(idx)}"

    def test_detects_reasonable_beat_count_360hz(self):
        ecg, expected = synthetic_ecg(360)
        amp, idx, delay = pan_tompkins(ecg, 360, gr=False)
        assert abs(len(idx) - expected) <= 2, f"Expected ~{expected} beats, got {len(idx)}"

    def test_detects_reasonable_beat_count_500hz(self):
        ecg, expected = synthetic_ecg(500)
        amp, idx, delay = pan_tompkins(ecg, 500, gr=False)
        assert abs(len(idx) - expected) <= 2, f"Expected ~{expected} beats, got {len(idx)}"


# ── Index validity ───────────────────────────────────────────────────────────

class TestIndexValidity:

    def test_indices_within_signal_bounds(self):
        fs = 200
        ecg, _ = synthetic_ecg(fs)
        amp, idx, delay = pan_tompkins(ecg, fs, gr=False)
        assert np.all(idx >= 0)
        assert np.all(idx < len(ecg))

    def test_indices_are_integers(self):
        ecg, _ = synthetic_ecg(200)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert idx.dtype in (np.int32, np.int64, int)

    def test_amplitudes_positive(self):
        ecg, _ = synthetic_ecg(200)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert np.all(amp >= 0)


# ── RR interval plausibility ─────────────────────────────────────────────────

class TestRRIntervals:

    def test_rr_intervals_plausible_200hz(self):
        """RR intervals should correspond to 40–200 bpm."""
        fs = 200
        ecg, _ = synthetic_ecg(fs, bpm=70)
        amp, idx, delay = pan_tompkins(ecg, fs, gr=False)
        if len(idx) > 1:
            rr = np.diff(idx) / fs
            assert np.all(rr > 0.3), "Some RR intervals imply >200 bpm"
            assert np.all(rr < 1.5), "Some RR intervals imply <40 bpm"

    def test_rr_intervals_plausible_360hz(self):
        fs = 360
        ecg, _ = synthetic_ecg(fs, bpm=70)
        amp, idx, delay = pan_tompkins(ecg, fs, gr=False)
        if len(idx) > 1:
            rr = np.diff(idx) / fs
            assert np.all(rr > 0.3)
            assert np.all(rr < 1.5)


# ── Delay ────────────────────────────────────────────────────────────────────

class TestDelay:

    def test_delay_non_negative(self):
        ecg, _ = synthetic_ecg(200)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert delay >= 0

    def test_delay_reasonable_magnitude(self):
        fs = 200
        ecg, _ = synthetic_ecg(fs)
        amp, idx, delay = pan_tompkins(ecg, fs, gr=False)
        assert delay <= fs


# ── Edge cases ───────────────────────────────────────────────────────────────

class TestEdgeCases:

    def test_input_accepted_as_list(self):
        ecg, _ = synthetic_ecg(200)
        amp, idx, delay = pan_tompkins(ecg.tolist(), 200, gr=False)
        assert len(idx) > 0

    def test_2d_input_raises(self):
        ecg = np.ones((2, 100))
        with pytest.raises((ValueError, Exception)):
            pan_tompkins(ecg, 200, gr=False)

    def test_high_bpm_signal(self):
        ecg, expected = synthetic_ecg(200, bpm=120)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert len(idx) > 0

    def test_low_bpm_signal(self):
        ecg, expected = synthetic_ecg(200, bpm=50)
        amp, idx, delay = pan_tompkins(ecg, 200, gr=False)
        assert len(idx) > 0