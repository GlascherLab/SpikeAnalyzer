# SpikeAnalyzer
MATLAB tool to detect unphysiological spikes in fMRI times series (uses SPM functions)

## RATIONALE 

Spike Analyzer is a program to help you to detect (unnatural) spike
(outliers) in an EPI time series. It is NOT a replacement for visual
inspection, but it'll help you in selecting certain suspicious scan for
closer visual inspection, e.g. with spm_movie.
Spikes are detected in two ways: (1) by thresholding the 1st derivative
of the global activity of an EPI series or (2) by thresholding the
movement velocity. All scans that exceed these thresholds are marked as
suspicious (candidates for removal by nuisance regressors). Images
detected by (1) are shown as vertical bars in light blue in the plots,
images detected by (2) are shown in light red.

Visually you can see spike as "jumps" in the global activity, or as
sudden jerks in the 1st derivative. The 1st derivative is better suited
for defining a threshold because except for some jerks, this curves will
be close to zero, even when the original global mean curve settles on a
different level after some jumps.

These "jumps" can be caused e.g. by scanner malfunction, but also by
subject movement. In order evaluate the latter cause, Spike Analyzer uses
the subject movement parameters to predict the global mean in a defined
windows of scans (defaults: 20). It then moves through the entire time
series with the same GLM. R-Square are color-coded and are used to
evaluate if some jumps in the global mean could have been caused by
movement (hot colorbar from dark red - orange - yellow - white). If the
jumps in the global mean curve coincide with a hot color in the R-Square
bar (named: Correlation between global activity and movement parameters)
then there is a good chance that it was cause by movement.
 
Movement parameters are also shown to evaluate, if the subject has shown
extraordinary large movement (rule of thumb: exclude subject if
tranlational movement exceed 2x voxel size or 1 degree of rotational
movement. These are independent of whether the movement is correlated
with the global mean.

The global mean is also shown with effects of movement removed (again
through a GLM). Jumps in the global mean curve that are correlated
with movement should not be visible in the last plot

Finally, Spike Analyzer will also shown you how many images are marked as
suspicious. You can then decide to adjust the thresholds (cfg.globthr and
cfg.movethr) in order to exclude or include more scans.

Spike Analyzer will output matrices of dummy regressors, in which a
suspicious scan is marked with 1 and all other scans are marked with
zero. The 4 matrices are for global mean based nuisance regressors
(gmNuis), movement velocity based nuisance regressors (mNuis), the
correlation-based neuisance regressors (rsqrNuis), or all combined
(allNuis). To include these nuisance regressors in your design 
matrix you need to have to enter them in SPM.Sess(i).C during design
specification, e.g. 

SPM.Sess(1).C.C = nuis.all;
SPM.Sess(1).C.name = repmat({'BS'},1,size(SPM.Sess(1).C.C,2));
(since one is usually not interested in looking at the effect of a single
nuisance regressor, they are all assigned the same name (BS = bad scan))

## EXPLANATION OF PLOTS (from top to bottom)

1. and 2. Normalized Global activity and 1st derivative of global
   activity (gray) 
   Configuration Options:
   cfg.globthr  = threshold (in SD) for excluding scans that exceed this
                  threshold from the 1st derivative of global mean
                  (default: 0.5) 
3. small bar with R-square values of predicting the normalized global
   mean from the movement parameters
   cfg.rsqrthr  = threshold (range 0-1) on R-Square value between
                  movement parameters as a predictor of global mean
                  (default: 0.8)

4. translational movement parameters
   cfg.ralignthr = threshold (1 or 2) for translational and rotational
   parameters (rule of thumb). Conservative: translational threshold
   1 * voxel size, rotational threshold: 1 degree. Liberal: 2* voxel
   size, 2 degrees
5. rotational movement parameters
6. Fast movement detection (using the algorithm of ArtRepair Toolbox)
   Configuration Options:
   cfg.movethr = threshold for excluding scan based on fast movement
                 (mm/TR) (default: 0.2)

7. Normalized Global Activity with effects of movement parameters removed
   (by multiple linear regression)

## SOME TECHNICAL INFORMATION

Spike Analyzer used to compute the global mean exclusively using 
spm_global. However, several test and simulations have revealed that 
spm_global is prone to include outlying voxels, which tend to distort the
global mean, thus leading to the potential exclusion of too many "bad
scans". spm_global computes the global mean on an image-by-image basis,
which almost certainly leads to an inclusion of different voxels for each
image. This is not optimal for the computation of a global mean time
series. Therefore, Spike Analyzer does not rely on spm_global anymore and
offers 3 options for selection the voxels to compute the global mean:
1. Segment the first EPI image of the series using spm_segment and
   computing a brain mask from the gray and white matter segments.
   Optionally, the gray and white matter segments can be used exclusively
   for voxel selection. This option is slow, but provides the most
   accurate voxel selection
2. Computing the global mean on a constant number of voxels. This
   requires several passes through the data: (a) compute and
   image-by-image threshold by dividing the mean of each image by 8 (like
   in the spm_global), (b) setting the overall threshold to the mean of
   the image-spcific thresholds, (c) finding the voxels in each image
   that exceed the overall threshold, and (d) selecting only those voxels
   common to all images that exceed the overall threshold. Although this
   seems a bit complicated, the results of this method are much better
   than the plain spm_global computations and are similar to the
   brainmask approach in 1.
3. Compute the global mean only in an anatomical ROI. If you are
   interested in a specific brain region, then you can configure an
   anatomical ROI image and only those voxels in the ROI will be selected
   for the computation of global mean activity.
4. For the sake of comparison, cfg.voxselect can be set to 'spm_global',
   which will then use the old spm_global approach. This method is
   depreciated.

The spm_global problem of selecting different voxels per image is
especially aggravated for images that are not preprocessed yet. It is
therefore recommended to use Spike Analyzer only on preprocessed images
that will enter the statistical analysis (i.e. realign and warped
images). The nuisance regressor created by Spike Analyzer will then fit
to the images that will be eneter into the 1st level analysis, which is
preferable. Happy Analyzing!
