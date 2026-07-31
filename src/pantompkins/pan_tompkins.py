"""
Pan-Tompkins QRS Detection Algorithm - Python Implementation

Complete implementation of the Pan-Tompkins algorithm for real-time QRS detection in ECG signals.

References:
    [1] Sedghamiz, H., "Matlab Implementation of Pan Tompkins ECG QRS detector.", 2014.
        https://www.researchgate.net/publication/313673153
    [2] Pan, J., Tompkins, W.J., "A Real-Time QRS Detection Algorithm",
        IEEE Transactions on Biomedical Engineering, VOL. BME-32, NO. 3, MARCH 1985.

Original MATLAB author: Hooman Sedghamiz, Feb 2018
                        MSc. Biomedical Engineering, Linköping University
                        Email: Hooman.sedghamiz@gmail.com

Python port:            Dr. Hatem Zehir
                        PhD. Biometrics, Badji Mokhtar - Annaba University
                        Email: hatem.zehir@gmail.com

Copyright (c) 2018, Hooman Sedghamiz
Copyright (c) 2026, Dr. Hatem Zehir
All rights reserved.

BSD 3-Clause License - see LICENSE file for full details.
"""

import numpy as np
from scipy.signal import butter, filtfilt, find_peaks
import matplotlib.pyplot as plt


def pan_tompkins(ecg: np.ndarray, fs: float, gr: bool = True):
    """
    Pan-Tompkins QRS Detection Algorithm.

    Parameters
    ----------
    ecg : np.ndarray
        Raw 1D ECG signal vector.
    fs : float
        Sampling frequency in Hz (e.g. 200, 360, 500).
    gr : bool, optional
        If True, plots intermediate and final results. Default is True.

    Returns
    -------
    qrs_amp_raw : np.ndarray
        Amplitudes of detected R waves (from bandpass-filtered signal).
    qrs_i_raw : np.ndarray
        Indices of detected R waves (from bandpass-filtered signal).
    delay : float
        Number of samples the signal is delayed due to filtering.
    """
    ecg = np.asarray(ecg, dtype=float).ravel()

    if ecg.ndim != 1:
        raise ValueError("ecg must be a 1D vector")

    # ── Initialize ──────────────────────────────────────────────────────────
    delay = 0
    skip = 0
    m_selected_RR = 0
    mean_RR = 0
    ser_back = 0

    # ── Noise cancellation / Bandpass filtering (5–15 Hz) ───────────────────
    if fs == 200:
        ecg = ecg - np.mean(ecg)

        # Low-pass ~12 Hz
        Wn_low = 12 * 2 / fs
        b_l, a_l = butter(3, Wn_low, btype='low')
        ecg_l = filtfilt(b_l, a_l, ecg)
        ecg_l = ecg_l / np.max(np.abs(ecg_l))

        # High-pass ~5 Hz
        Wn_high = 5 * 2 / fs
        b_h, a_h = butter(3, Wn_high, btype='high')
        ecg_h = filtfilt(b_h, a_h, ecg_l)
        ecg_h = ecg_h / np.max(np.abs(ecg_h))

        if gr:
            fig1, axes = plt.subplots(3, 2, figsize=(12, 8))
            axes[0, 0].plot(ecg); axes[0, 0].set_title('Raw Signal'); axes[0, 0].axis('tight')
            axes[0, 1].plot(ecg_l); axes[0, 1].set_title('Low Pass Filtered'); axes[0, 1].axis('tight')
            axes[1, 0].plot(ecg_h); axes[1, 0].set_title('High Pass Filtered'); axes[1, 0].axis('tight')
    else:
        f1, f2 = 5, 15
        Wn = [f1 * 2 / fs, f2 * 2 / fs]
        b_bp, a_bp = butter(3, Wn, btype='band')
        ecg_h = filtfilt(b_bp, a_bp, ecg)
        ecg_h = ecg_h / np.max(np.abs(ecg_h))

        if gr:
            fig1, axes = plt.subplots(3, 2, figsize=(12, 8))
            axes[0, 0].plot(ecg); axes[0, 0].set_title('Raw Signal'); axes[0, 0].axis('tight')
            axes[1, 0].plot(ecg_h); axes[1, 0].set_title('Band Pass Filtered'); axes[1, 0].axis('tight')

    # ── Derivative filter ────────────────────────────────────────────────────
    # H(z) = (1/8T)(−z^−2 − 2z^−1 + 2z + z^2)
    if fs != 200:
        int_c = (5 - 1) / (fs * 1 / 40)
        b_d = np.interp(
            np.arange(1, 5 + 1e-10, int_c),
            np.arange(1, 6),
            np.array([1, 2, 0, -2, -1]) * (1 / 8) * fs
        )
    else:
        b_d = np.array([1, 2, 0, -2, -1]) * (1 / 8) * fs

    ecg_d = filtfilt(b_d, [1], ecg_h)
    ecg_d = ecg_d / np.max(np.abs(ecg_d))

    if gr:
        axes[1, 1].plot(ecg_d); axes[1, 1].set_title('Derivative Filtered'); axes[1, 1].axis('tight')

    # ── Squaring ─────────────────────────────────────────────────────────────
    ecg_s = ecg_d ** 2

    if gr:
        axes[2, 0].plot(ecg_s); axes[2, 0].set_title('Squared'); axes[2, 0].axis('tight')

    # ── Moving average (150 ms window) ───────────────────────────────────────
    win = round(0.150 * fs)
    ecg_m = np.convolve(ecg_s, np.ones(win) / win, mode='full')
    delay += win / 2

    if gr:
        axes[2, 1].plot(ecg_m)
        axes[2, 1].set_title('Moving Average (~150 ms window)')
        axes[2, 1].axis('tight')
        plt.tight_layout()
        plt.show()

    # ── Fiducial marks ───────────────────────────────────────────────────────
    min_dist = round(0.2 * fs)
    locs, props = find_peaks(ecg_m, distance=min_dist)
    pks = ecg_m[locs]
    LLp = len(pks)

    # ── Buffer initialisation ────────────────────────────────────────────────
    qrs_c = np.zeros(LLp)
    qrs_i = np.zeros(LLp, dtype=int)
    qrs_i_raw = np.zeros(LLp, dtype=int)
    qrs_amp_raw = np.zeros(LLp)
    nois_c = np.zeros(LLp)
    nois_i = np.zeros(LLp, dtype=int)
    SIGL_buf = np.zeros(LLp)
    NOISL_buf = np.zeros(LLp)
    THRS_buf = np.zeros(LLp)
    SIGL_buf1 = np.zeros(LLp)
    NOISL_buf1 = np.zeros(LLp)
    THRS_buf1 = np.zeros(LLp)

    # ── Training phase (first 2 s) ───────────────────────────────────────────
    train = ecg_m[:2 * int(fs)]
    THR_SIG = np.max(train) / 3
    THR_NOISE = np.mean(train) / 2
    SIG_LEV = THR_SIG
    NOISE_LEV = THR_NOISE

    train_h = ecg_h[:2 * int(fs)]
    THR_SIG1 = np.max(train_h) / 3
    THR_NOISE1 = np.mean(train_h) / 2
    SIG_LEV1 = THR_SIG1
    NOISE_LEV1 = THR_NOISE1

    Beat_C = 0
    Beat_C1 = 0
    Noise_Count = 0

    # ── Main decision loop ───────────────────────────────────────────────────
    for i in range(LLp):
        # Locate corresponding peak in bandpass-filtered signal
        start_idx = locs[i] - round(0.150 * fs)
        if start_idx >= 1 and locs[i] <= len(ecg_h):
            segment = ecg_h[start_idx: locs[i]]
            x_i = int(np.argmax(segment))
            y_i = segment[x_i]
        else:
            if i == 0:
                segment = ecg_h[:locs[i]]
                x_i = int(np.argmax(segment))
                y_i = segment[x_i]
                ser_back = 1
            elif locs[i] >= len(ecg_h):
                start_idx = max(0, locs[i] - round(0.150 * fs))
                segment = ecg_h[start_idx:]
                x_i = int(np.argmax(segment))
                y_i = segment[x_i]
            else:
                y_i, x_i = 0.0, 0

        # Update mean RR interval after 9 beats
        if Beat_C >= 9:
            diffRR = np.diff(qrs_i[Beat_C - 8: Beat_C])
            mean_RR = np.mean(diffRR)
            comp = qrs_i[Beat_C - 1] - qrs_i[Beat_C - 2]

            if comp <= 0.92 * mean_RR or comp >= 1.16 * mean_RR:
                m_selected_RR = m_selected_RR  # irregular → don't update
                THR_SIG = 0.5 * THR_SIG
                THR_SIG1 = 0.5 * THR_SIG1
            else:
                m_selected_RR = mean_RR

        # Choose test interval
        if m_selected_RR:
            test_m = m_selected_RR
        elif mean_RR and m_selected_RR == 0:
            test_m = mean_RR
        else:
            test_m = 0

        not_nois = 0

        # Search-back for missed QRS
        if test_m and Beat_C >= 1:
            if (locs[i] - qrs_i[Beat_C - 1]) >= round(1.66 * test_m):
                search_start = qrs_i[Beat_C - 1] + round(0.200 * fs)
                search_end = locs[i] - round(0.200 * fs)
                if search_start < search_end:
                    seg = ecg_m[search_start: search_end]
                    local_max_idx = int(np.argmax(seg))
                    pks_temp = seg[local_max_idx]
                    locs_temp = search_start + local_max_idx

                    if pks_temp > THR_NOISE:
                        qrs_c[Beat_C] = pks_temp
                        qrs_i[Beat_C] = locs_temp
                        Beat_C += 1

                        end_bp = min(locs_temp, len(ecg_h))
                        start_bp = max(0, locs_temp - round(0.150 * fs))
                        seg_bp = ecg_h[start_bp: end_bp]
                        if len(seg_bp):
                            x_i_t = int(np.argmax(seg_bp))
                            y_i_t = seg_bp[x_i_t]
                        else:
                            y_i_t, x_i_t = 0.0, 0

                        if y_i_t > THR_NOISE1:
                            qrs_i_raw[Beat_C1] = start_bp + x_i_t
                            qrs_amp_raw[Beat_C1] = y_i_t
                            Beat_C1 += 1
                            SIG_LEV1 = 0.25 * y_i_t + 0.75 * SIG_LEV1

                        not_nois = 1
                        SIG_LEV = 0.25 * pks_temp + 0.75 * SIG_LEV

        # Classify peak
        if pks[i] >= THR_SIG:
            # T-wave check
            if Beat_C >= 3:
                if (locs[i] - qrs_i[Beat_C - 1]) <= round(0.360 * fs):
                    s1 = locs[i] - round(0.075 * fs)
                    s2 = qrs_i[Beat_C - 1] - round(0.075 * fs)
                    Slope1 = np.mean(np.diff(ecg_m[max(0, s1): locs[i]]))
                    Slope2 = np.mean(np.diff(ecg_m[max(0, s2): qrs_i[Beat_C - 1]]))
                    if abs(Slope1) <= abs(0.5 * Slope2):
                        nois_c[Noise_Count] = pks[i]
                        nois_i[Noise_Count] = locs[i]
                        Noise_Count += 1
                        skip = 1
                        NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1
                        NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV
                    else:
                        skip = 0

            if skip == 0:
                qrs_c[Beat_C] = pks[i]
                qrs_i[Beat_C] = locs[i]
                Beat_C += 1

                if y_i >= THR_SIG1:
                    if ser_back:
                        qrs_i_raw[Beat_C1] = x_i
                    else:
                        qrs_i_raw[Beat_C1] = locs[i] - round(0.150 * fs) + x_i
                    qrs_amp_raw[Beat_C1] = y_i
                    Beat_C1 += 1
                    SIG_LEV1 = 0.125 * y_i + 0.875 * SIG_LEV1

                SIG_LEV = 0.125 * pks[i] + 0.875 * SIG_LEV

        elif THR_NOISE <= pks[i] < THR_SIG:
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV

        elif pks[i] < THR_NOISE:
            nois_c[Noise_Count] = pks[i]
            nois_i[Noise_Count] = locs[i]
            Noise_Count += 1
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV

        # Adjust thresholds
        if NOISE_LEV != 0 or SIG_LEV != 0:
            THR_SIG = NOISE_LEV + 0.25 * abs(SIG_LEV - NOISE_LEV)
            THR_NOISE = 0.5 * THR_SIG

        if NOISE_LEV1 != 0 or SIG_LEV1 != 0:
            THR_SIG1 = NOISE_LEV1 + 0.25 * abs(SIG_LEV1 - NOISE_LEV1)
            THR_NOISE1 = 0.5 * THR_SIG1

        SIGL_buf[i] = SIG_LEV
        NOISL_buf[i] = NOISE_LEV
        THRS_buf[i] = THR_SIG
        SIGL_buf1[i] = SIG_LEV1
        NOISL_buf1[i] = NOISE_LEV1
        THRS_buf1[i] = THR_SIG1

        skip = 0
        not_nois = 0
        ser_back = 0

    # ── Trim output arrays ───────────────────────────────────────────────────
    qrs_i_raw = qrs_i_raw[:Beat_C1]
    qrs_amp_raw = qrs_amp_raw[:Beat_C1]
    qrs_c = qrs_c[:Beat_C]
    qrs_i = qrs_i[:Beat_C]

    # ── Plotting ─────────────────────────────────────────────────────────────
    if gr:
        # Overlay on moving-average signal
        plt.figure(figsize=(12, 4))
        plt.plot(ecg_m, label='MVI Signal')
        plt.scatter(qrs_i, qrs_c, color='m', zorder=5, label='QRS')
        plt.plot(locs, NOISL_buf, '--k', linewidth=2, label='Noise Level')
        plt.plot(locs, SIGL_buf, '--r', linewidth=2, label='Signal Level')
        plt.plot(locs, THRS_buf, '--g', linewidth=2, label='Threshold')
        plt.title('QRS on MVI Signal')
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Summary figure
        fig2, az = plt.subplots(3, 1, figsize=(12, 9))
        az[0].plot(ecg_h); az[0].set_title('QRS on Filtered Signal'); az[0].axis('tight')
        az[0].scatter(qrs_i_raw, qrs_amp_raw, color='m', zorder=5)
        az[0].plot(locs, NOISL_buf1, '--k', linewidth=2)
        az[0].plot(locs, SIGL_buf1, '-.r', linewidth=2)
        az[0].plot(locs, THRS_buf1, '-.g', linewidth=2)

        az[1].plot(ecg_m); az[1].set_title('QRS on MVI Signal'); az[1].axis('tight')
        az[1].scatter(qrs_i, qrs_c, color='m', zorder=5)
        az[1].plot(locs, NOISL_buf, '--k', linewidth=2)
        az[1].plot(locs, SIGL_buf, '-.r', linewidth=2)
        az[1].plot(locs, THRS_buf, '-.g', linewidth=2)

        ecg_zero = ecg - np.mean(ecg)
        az[2].plot(ecg_zero); az[2].set_title('QRS Pulse Train on Raw ECG'); az[2].axis('tight')
        for idx in qrs_i_raw:
            az[2].axvline(x=idx, color='r', linewidth=1.5, linestyle='-.')
        
        plt.tight_layout()
        plt.show()

    return qrs_amp_raw, qrs_i_raw, delay


# ── Quick smoke-test ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Generate a simple synthetic ECG-like signal for testing
    fs = 200
    t = np.linspace(0, 10, int(10 * fs))

    # Approximate QRS complexes as narrow Gaussians at ~1 Hz
    ecg_test = np.zeros_like(t)
    beat_times = np.arange(0.5, 10, 0.85)  # ~70 bpm
    for bt in beat_times:
        ecg_test += np.exp(-((t - bt) ** 2) / (2 * (0.02) ** 2))

    # Add some noise
    rng = np.random.default_rng(42)
    ecg_test += 0.05 * rng.standard_normal(len(t))

    amp, idx, delay = pan_tompkins(ecg_test, fs, gr=True)

    print(f"Detected {len(idx)} QRS complexes")
    print(f"R-peak indices: {idx}")
    print(f"Delay: {delay:.1f} samples")