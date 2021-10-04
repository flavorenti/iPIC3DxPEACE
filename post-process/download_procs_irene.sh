#!/bin/bash

## PR0        --> run1 --> 5500
## PR1        --> run1 --> 7500
## PR2        --> run28--> 11200
## PR0-bigbox --> run1 --> 9500


while read line; do 
echo $line
rsync --progress -av lavorenf@irene-eu.ccc.cea.fr:/ccc/scratch/cont005/gen12622/lavorenf/Mercury_SaeInit/PR1/run1/data/restart$line.hdf ~/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/run1/data/ 
done < good_patches_MarinerX_PR1.txt


