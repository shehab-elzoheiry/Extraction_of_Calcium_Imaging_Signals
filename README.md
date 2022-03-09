# Extraction_of_Calcium_Imaging_Signals

## Content:

This repository includes the functions used to extract calcium imaging signals from several regions of interest. Main purpose is to detect "clean" signal after discarding contamination due to cross-talk, which is common for wide field / fluorescent microscopes.

This work is part of the publication [Elzoheiry et al 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7820691/pdf/10.1177_0271678X19892657.pdf). Check the methods part and supplementary Figure.3 for more details. Below is brief explanation.

## Explanation of decontamination:

For a true signal to be detected (Supplementary Figure 3(a)), two conditions have to be satisfied. **First**, traces from all pixels within the relevant ROI should show higher values (at least 95% of the pixels) than the baseline fit of the traces. This indicates that the whole ROI is fluorescing which occurs if a cell is active. However, out of focus background can increase pixel values in the ROI as well, leading to false positive signal. Therefore, the **second** condition is that the average pixel values from the surrounding region (border of the ROI) are lower than the values of pixels within the ROI (less than the lowest 5%).

On the contrary, when contamination of the complete ROI occurs (Supplementary Figure 3(b)), transients from the surrounding region are high enough relative to the values of pixels within the ROI (higher than the lowest 5%), thereby, violating the second condition described above. In this case, the difference between activity outside and within the ROI during the events was linearly interpolated between the beginning and end of the relevant event. In the case of partial overlap contamination (Supplementary Figure 3(c)), violation of the first condition occurs, in which only a portion of pixels show elevated transients. These pixels are clustered using K-means and their fluorescent transients are discarded from calculating the final average fluorescent signal.

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

## Implementation on local machine:
* Before running extract_denoise_calcium_imaging.m, you should generate ROIs (regions of interest) using imageJ or fiji.
* Determine the directories for the .tif files, ROIs, and the directory where you want to save the output
* run the extract_denoise_calcium_imaging.m function.
* the remaining two functions are generating:
    a. Border surrounding each ROI
    b. Baseline fit
  Both are required for deduction of background and cross-talk contamination
