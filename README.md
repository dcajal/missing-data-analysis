# missing-data-analysis

MATLAB code used in Cajal, D., Hernando, D., Lázaro, J., Laguna, P., Gil, E., & Bailón, R. (2022). Effects of missing data on heart rate variability metrics. Sensors, 22(15), 5774.
See: https://www.mdpi.com/1424-8220/22/15/5774

Heart rate variability (HRV) has been studied for decades in clinical environments. Currently, the exponential growth of wearable devices in health monitoring is leading to new challenges that need to be solved. These devices have relatively poor signal quality and are affected by numerous motion artifacts, with data loss being the main stumbling block for their use in HRV analysis. In the present paper, it is shown how data loss affects HRV metrics in the time domain and frequency domain and Poincaré plots. A gap-filling method is proposed and compared to other existing approaches to alleviate these effects, both with simulated (16 subjects) and real (20 subjects) missing data. Two different data loss scenarios have been simulated: (i) scattered missing beats, related to a low signal to noise ratio; and (ii) bursts of missing beats, with the most common due to motion artifacts. In addition, a real database of photoplethysmography-derived pulse detection series provided by Apple Watch during a protocol including relax and stress stages is analyzed. The best correction method and maximum acceptable missing beats are given. Results suggest that correction without gap filling is the best option for the standard deviation of the normal-to-normal intervals (SDNN), root mean square of successive differences (RMSSD) and Poincaré plot metrics in datasets with bursts of missing beats predominance (p<0.05), whereas they benefit from gap-filling approaches in the case of scattered missing beats (p<0.05). Gap-filling approaches are also the best for frequency-domain metrics (p<0.05). The findings of this work are useful for the design of robust HRV applications depending on missing data tolerance and the desired HRV metrics.

Database not included.
