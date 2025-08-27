# Stretching-optimization

[English](README.md) | [中文](README_zh.md)

***

As an English learner, I still have a poor command of English. Thus, this document may inevitably contain some inaccuracies. Hopefully, that will not mislead you.

***

[MSNoise1.6](http://www.msnoise.org ) "A Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise"

This document is on the optimization of the stretching code in [MSNoise1.6](http://www.msnoise.org ), particularly focusing on improving the ST (Stretching) algorithm by incorporating clock error correction. This is just a simple attempt, and the code remains rudimentary, with ideas yet to be fully validated. Fortunately, the code may not report error, if the paths are correctly configured.

### Debug

2025 04 19
The original ST-FPC methodology employed a folding-based approach for MSE (Mean Squared Error) computation, which introduced structural inconsistencies and led to significant information degradation. In contrast, the upgraded ST-PC algorithm yields smoother operational profiles while fixing a bug in the legacy elimination module. Nevertheless, I maintain the recommendation to prioritize the updated elimination method for outlier rejection tasks.

<img src="Figure/smooth.jpg?v=1&type=image" width="500">

Btw a 20DAYS stack with ST_FPC suggests great advantages

<img src="Figure/OFFSET.png?v=1&type=image" width="500">

### Usage

The new ST_PC(Prediction-Correction）is highly recommended for archiving smoother velocity transitions and enhanced mathematical validity, sharing the same usage with the ST_FPC.

The ST_FPC (Fold-Prediction-Correction) is designed for the stretching calculation step after stacking in the MSNoise workflow. It integrates with MSNoise by using the STACKS directory as input. Parameters such as file path, lag time window, stretching range, and step size must be reconfigured in the code header.
  
* readsac.m is an official SAC function, which should share the same folder with ST_FPC.m .

* stretch.py is just the MSNoise for ST calculations.

* MATLAB requires a Python environment to call the map_coordinates function.

### Code Explanation

In my view, the MWCS (Moving Window Cross-Spectral) method is theoretically superior to ST (Stretching) for calculating velocity changes ( $dv/v=-dt/t$) between reference (Ref) and daily (Days) waveforms. MWCS directly utilizes phase relationships to fit dt/t, whereas ST relies on the correlation coefficient (CC, Pearson) to determine dt/t(Obermann and Hillers, 2019), introducing two layers of distortion: from waveform to CC, and from CC to dt/t. However, ST involves fewer parameters and is simpler to implement, making it a practical reference.

The calculation of velocity changes $dt/t$ from empirical Green’s functions involves determining the stretching coefficient $a$ between waveforms Ref and Days, modeled as $f(x)$ vs. $f(ax + b)$. Clock error manifests as the offset $b$ in this framework. Below, we discuss dt/t and offset separately, using test data from six nearby stations over several months.

### dt/t

First, the Ref waveform is stretched over a predefined range with incremental steps. For each stretched Ref waveform, the correlation coefficient with the Days waveform is computed. The stretching coefficient corresponding to the maximum CC is selected as dt/t.

Figure 1 shows a cross-correlation waveform (-120–120 s), with the yellow region marking the lag time window (20–80 s). During stretching, the original code(MSNoise1.6 stretching.py) selects data points within the yellow window (including non-window regions with zeroing), stretches the time axis via scipy.ndimage.map_coordinates, and computes CC using all stretched data points (-120–120 s) against the Days waveform (non-window regions zeroed).

<img src="Figure/c755c4a0-e22c-11ef-b911-b3d360b4824a.jpeg?v=1&type=image" width="500">

In the first place, replicate the original Python code with MATLAB. However, MATLAB’s interp1 lacks the pre-filtering step comparing with map_coordinates, leading to inferior results. To address this, the Python function was directly called. The results, Figure 2 (Colors represent the magnitude of CC) closely match the original Python implementation.

<img src="Figure/8dbdb1c0-e2d7-11ef-bda0-e3f385aefa20.jpeg?v=1&type=image" width="500">

It seems that there are some issues within the original code:

When compressing waveforms, edge values are defaulted as zeros. However, the full waveform are available, so fixed window can be extended to ensure valid interpolation.

When stretching waveforms, additional data points outside the window are generated but not removed. In addition, there is a zero interval between the two data segments, and the corresponding test part is [-20, 20]. These may all lead to impact CC calculations.

For the waveform segment, I prefer to extend proper window, and then fold the double-sides waveform as the computeds waveform segment

Figure 3 shows that selecting valid segments increased CC values by ~10% for 89% of cases.

<img src="Figure/31db0ce0-e2d7-11ef-bda0-e3f385aefa20.jpeg?v=1&type=image" width="500">

Figure 4 demonstrates that window extension slightly improved CC for 85% of cases.

<img src="Figure/0ce2d0f0-e235-11ef-b911-b3d360b4824a.jpeg?v=1&type=image" width="500">

### Offset

Unlike MWCS(Stehly et al., 2007), ST lacks a robust method to handle clock error. A simple solution is to fit Ref and Days waveforms to some form of function expression, like $f(x)$ vs. $f(ax + b)$, and then compare their coefficients for $b$. However, directly aiming at the targret window seems difficult. After some attempts, this Gaussian fit based prediction-correction method has been developed.

Figure 5 shows the cross-correlation functions (yellow: Days, purple: Ref) for a station pair. As the full waveform resembles Gaussian distribution, taking the absolute value of the data and the difference of two centers (derived from 4th-order Gaussian fitting) suggests an estimated offset.

<img src="Figure/01d53ce0-e2d8-11ef-bda0-e3f385aefa20.jpeg?v=1&type=image" width="500">

Figure 6 illustrates the horizontal distribution of offset derived from 4th-order Gaussian fitting across 15 station pairs, indicating correlations with station pairs. This kind offset may be dependable as the alternative of clock error.

<img src="Figure/9a4db2d0-e473-11ef-ad60-9304b2563e4d.jpeg?v=1&type=image" width="500">

Figure 7 There are 15 pairs of stations on the vertical axis, and horizontal bands with similar colors can be seen, indicating that the offset results are somewhat correlated with the pairs of stations

<img src="Figure/0b0b7de0-e46f-11ef-ad60-9304b2563e4d.jpeg?v=1&type=image" width="500">

It should be noted that the current offset has no direct correlation with waveform phase relationship, being merely an artifact of the simulation. Moreover, Figure 8, the obtained CC results, demonstrate only marginal improvement in approximately half of the cases, exhibiting insufficient reliability. Therefore, the selection of optimal offsets must be integrated with CC methodology.

<img src="Figure/1f9fcc70-e52d-11ef-baf1-a700618c26d1.jpeg?v=1&type=image" width="500">

Additionally, proper fitting parameters are crucial to Gaussian fitting. Figure 9 demonstrates the impact of poor Gaussian fitting parameters on offset accuracy.

<img src="Figure/438147e0-e474-11ef-ad60-9304b2563e4d.jpeg?v=1&type=image" width="500">

However, if just depend on maxCC, the result of offset will show the characteristic of random distribution, like Figure 10

<img src="Figure/1c80f5d0-e886-11ef-8864-3d63cc015e26.jpeg?v=1&type=image" width="500">

A re-examination of its correlation with station distribution revealed extensive color interlacing zones, Figure 11, indicating oscillatory regions with significant data fluctuations. The offset in these segments cannot be attributed to clock errors.

<img src="Figure/83dde800-e52a-11ef-baf1-a700618c26d1.jpeg?v=1&type=image" width="500">

Figure 12 shows improved offset after oscillation removal and boundary constraints.

<img src="Figure/96fe9370-e887-11ef-8864-3d63cc015e26.jpeg?v=1&type=image" width="500">

#### Considering above issues, the prediction-correction is recommended. Figure 13 illustrates the full workflow.

#### 1. Prediction: Gaussian fitting predicts an initial pre-offset. Comopute pre-dt/t with pre-offset.

#### 2. Correction: Depending on the stretched waveform, corr-offset is determined by maximizing CC within the predicted offset's range. Eliminate the oscillation zone.

#### 3. With final corr-offset, recompute the corr-dt/t.

<img src="Figure/624f6590-eb7b-11ef-bb64-49bc940d885a.jpeg?v=1&type=image" width="500">

Figure 14 demonstrates that clock error correction increased CC by ~15% for 88% of cases compared to the Fold method.

<img src="Figure/5ae41860-e2c4-11ef-b8f7-7fbcf1303d2b.jpeg?v=1&type=image" width="500">

Finnally, Figure 15 compares CC improvements: 97% of cases show ~15% average increase.

<img src="Figure/cf4b1880-e2c3-11ef-b8f7-7fbcf1303d2b.jpeg?v=1&type=image" width="500">

Based on imporved CC, Figure 16 shows stable velocity changes if using a 3-point weighted average of CC.

<img src="Figure/914b5820-eba5-11ef-b0c0-5937d7efb862.jpeg?v=1&type=image" width="500">

### Conclusion

Refer to MSNoise1.6 stretching, by modifying the waveform segment (Fold), Gaussian fitting, prediction-correction, and oscillation filtering, yield higher CC values for dt/t. However, challenges remain:

1. How does ST’s offset compare to MWCS’s clock error correction?

2. How to optimize Gaussian fitting parameters or identify better fitting methods?

3. Are there better approaches than zeroing oscillatory regions?

4. The offset distributions appear to kind of longitudinal correlation, or say periodic pattern. Can this periodic pattern be modeled for further correction?

5. Would alternative correlation metrics (e.g., distance correlation, MIC) outperform Pearson’s CC?

#####  My personal study on ambient noise is still relatively shallow, and my understanding of related issues is not yet deep. These questions warrant further investigation. For any error, including issues in language expression, present in this document, I look forward to your corrections!

***

References

[1]Obermann A, Hillers G. Chapter Two-Seismic time-lapse interferometry across scales [J]. Advance in Geophysics, 2019, 60: 65-143.

[2]Stehly L., Campillo M., Shapiro N. M. 2007. Traveltime measurements from noise correlation: stability and detection of instrumental time-shifts. Geophys. J. Int., 171(4): 223-230.
