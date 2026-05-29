# Pan-Tompkins QRS Detector — Python

A Python port of the Pan-Tompkins real-time QRS detection algorithm for ECG signals.

> **Original MATLAB implementation** by Hooman Sedghamiz (Feb 2018), MSc Biomedical Engineering, Linköping University.  
> **Python port** by Dr. Hatem Zehir.  
> Python conversion retains the original BSD 3-Clause license — see [`LICENSE`](LICENSE).

---

## Algorithm Overview

The Pan-Tompkins algorithm is a classic, widely used method for detecting QRS complexes (R-peaks) in ECG signals in real time. It processes the ECG through a signal-processing pipeline:

1. **Bandpass filter** (5–15 Hz) — removes baseline wander and high-frequency noise
2. **Derivative filter** — emphasises the steep slopes of QRS complexes
3. **Squaring** — makes all values positive and amplifies large slopes
4. **Moving-window integration** (~150 ms window) — produces a smooth envelope
5. **Adaptive thresholding & decision logic** — distinguishes true QRS complexes from T-waves and noise using two running thresholds updated after each beat

### References
- Pan, J. & Tompkins, W. J., *"A Real-Time QRS Detection Algorithm"*, IEEE Trans. Biomed. Eng., BME-32(3), March 1985.
- Sedghamiz, H., *"Matlab Implementation of Pan Tompkins ECG QRS detector"*, 2014. [ResearchGate](https://www.researchgate.net/publication/313673153)

---

## Installation

```bash
pip install -r requirements.txt
```

**Dependencies:**

| Package | Version |
|---------|---------|
| numpy | ≥ 1.21 |
| scipy | ≥ 1.7 |
| matplotlib | ≥ 3.4 |

---

## Usage

```python
import numpy as np
from pan_tompkins import pan_tompkins

# Load your ECG signal and sampling frequency
# ecg : 1-D numpy array of the raw ECG
# fs  : sampling frequency in Hz (e.g. 200, 360, 500)

qrs_amp, qrs_idx, delay = pan_tompkins(ecg, fs, gr=True)

print(f"Detected {len(qrs_idx)} QRS complexes")
print(f"R-peak sample indices: {qrs_idx}")
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `ecg` | `np.ndarray` | Raw 1-D ECG signal |
| `fs` | `float` | Sampling frequency in Hz |
| `gr` | `bool` | Plot intermediate stages and results (default `True`) |

### Returns

| Return value | Description |
|---|---|
| `qrs_amp_raw` | Amplitudes of detected R-waves (from bandpass-filtered signal) |
| `qrs_i_raw` | Sample indices of detected R-waves |
| `delay` | Filter delay in samples |

---

## Quick Demo

Running the script directly generates a synthetic ECG (~70 bpm, 10 s at 200 Hz) and detects QRS complexes:

```bash
python pan_tompkins.py
```

Expected output:
```
Detected 12 QRS complexes
R-peak indices: [ 100  270  440  610  780  950 1119 1289 1459 1629 1799 1969]
Delay: 15.0 samples
```

---

## Notes

- For `fs = 200 Hz`, separate low-pass (12 Hz) and high-pass (5 Hz) Butterworth filters are applied sequentially.
- For all other sampling rates, a single bandpass Butterworth filter (5–15 Hz) is used.
- The derivative filter kernel is interpolated to match the sampling rate when `fs ≠ 200`.
- `filtfilt` (zero-phase filtering) is used throughout, so `delay` reports the moving-average window contribution only.

---

## License

BSD 3-Clause License. Copyright (c) 2018, Hooman Sedghamiz. See [`LICENSE`](LICENSE) for full details.

Copyright (c) 2018 Hooman Sedghamiz  
Copyright (c) 2026 Dr. Hatem Zehir

---

## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome!

Please open an issue or submit a pull request.
