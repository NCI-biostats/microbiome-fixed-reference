# Microbiome-Fixed-Reference

## Description

Using standard microbiome reference groups to simplify beta-diversity analyses and faclitate independent validation.

## Functions

* code_microbiome_helper_function.r:  All the functions needed for running analyses on microbiome data.
* application/hmp_mean_stool_nasal.r:  For HMP project, calculate distance between reference and test set. Reference sets are nasal reference and stool reference.
* application/hmp_median_stool_nasal.r: For HMP project, calculate median distance between reference and test set. Reference sets are stool reference and nasal reference.
* application/half_hmp_median_mean_stool_nasal.r: For HMP project, using 25% samples as reference set, 75% samples as test set. Calculate mean and median distance vs stool and nasal reference.
* application/hmp_ag_combine.r: Provide combined data sets for both HMP project and AG project.
* application/hmp_distance_dataset.r: Provide distance w.r.t to HMP stool reference and nasal reference.



