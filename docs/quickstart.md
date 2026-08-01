# Quickstart Tutorial

This guide shows you how to install `pantompkins` and run your first QRS detection in under five minutes.

---

## Installation

Install from PyPI:

```bash
pip install pantompkins
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/Hatem-Zehir/pan-tompkins-qrs-detector.git
```

---

## Basic Usage

```python
import numpy as np
from pantompkins import pan_tompkins

# Load your ECG signal (1-D NumPy array) and sampling frequency
ecg = np.loadtxt("my_ecg.csv")   # replace with your data
fs = 360                          # Hz — adjust to match your recording

# Detect QRS complexes
qrs_amp, qrs_idx, delay = pan_tompkins(ecg, fs, gr=True)

print(f"Detected {len(qrs_idx)} QRS complexes")
print(f"R-peak indices (samples): {qrs_idx}")
print(f"R-peak times (seconds):   {qrs_idx / fs}")
```

Setting `gr=True` (default) displays six figures showing each stage of the processing
pipeline. Set `gr=False` to suppress all plots.

---

## Working with Real ECG Data (PhysioNet MIT-BIH)

The [MIT-BIH Arrhythmia Database](https://physionet.org/content/mitdb/1.0.0/) is a free,
widely used benchmark. Load it with the `wfdb` library:

```bash
pip install wfdb
```

```python
import wfdb
from pantompkins import pan_tompkins

# Stream record 100 from MIT-BIH (no local download needed)
record = wfdb.rdrecord("100", sampfrom=0, sampto=10800, pn_dir="mitdb")

ecg = record.p_signal[:, 0]   # Lead MLII
fs  = record.fs                # 360 Hz

qrs_amp, qrs_idx, delay = pan_tompkins(ecg, fs, gr=True)
print(f"Detected {len(qrs_idx)} beats in 30 s")
```

---

## Computing Heart Rate

```python
rr_intervals   = np.diff(qrs_idx) / fs   # seconds
heart_rate_bpm = 60 / rr_intervals
print(f"Mean HR: {heart_rate_bpm.mean():.1f} bpm")
```

---

## Function Reference

### `pan_tompkins(ecg, fs, gr=True)`

**Parameters**

| Name  | Type          | Description                               |
|-------|---------------|-------------------------------------------|
| `ecg` | `array-like`  | Raw 1-D ECG signal                        |
| `fs`  | `float`       | Sampling frequency in Hz                  |
| `gr`  | `bool`        | Plot intermediate stages (default `True`) |

**Returns**

| Name           | Type          | Description                                      |
|----------------|---------------|--------------------------------------------------|
| `qrs_amp_raw`  | `np.ndarray`  | R-wave amplitudes (bandpass-filtered signal)     |
| `qrs_i_raw`    | `np.ndarray`  | R-wave sample indices                            |
| `delay`        | `float`       | Filter delay in samples                          |

---

## Processing Pipeline

The algorithm applies five sequential stages:

1. **Bandpass filter (5–15 Hz)** — removes baseline wander and high-frequency noise
2. **Derivative filter** — emphasises steep QRS slopes
3. **Squaring** — amplifies large slopes, suppresses P/T waves
4. **Moving-window integration (~150 ms)** — smooths the signal envelope
5. **Adaptive thresholding** — classifies peaks with T-wave rejection and search-back

See the [original paper](https://doi.org/10.1109/TBME.1985.325532) for full details.