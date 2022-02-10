#!/bin/bash
##########################################################
############### START POST-PROCESS #######################
##########################################################
export NAME=test_simu
export WORK=$WORK
export SCRATCH=$SCRATCH
declare -i nrun=1
##########################################################
##########################################################

#-- MAKE DIRECTORIES
mkdir $WORK/$NAME
mkdir $WORK/$NAME/images
mkdir $WORK/$NAME/texts

#-- COMPUTATIONAL TIME
python print_iPIC3D_comp-time.py $SCRATCH/$NAME $WORK/$NAME/texts/
python plot_comp-time_iPIC3D.py $WORK/$NAME

#-- ENERGY EVOLUTION
for irun in $(seq 0 $nrun)
do
    cp $SCRATCH/$NAME/run$irun/data/ConservedQuantities.txt $WORK/$NAME/texts/ConservedQuantities$irun.txt
done
python plot_iPIC3D_energy.py $SCRATCH/$NAME $WORK/$NAME/images

#-- MEMORY USAGE
echo "# ARRAY FIELDS MEMORY [KB]:" > $WORK/$NAME/texts/memory.txt
du -h $SCRATCH/$NAME/run0/data/*_B_0.vtk >> $WORK/$NAME/texts/memory.txt
du    $SCRATCH/$NAME/run0/data/*_B_0.vtk >> $WORK/$NAME/texts/memory.txt
echo "# SCALAR FIELD MEMORY [KB]:" >> $WORK/$NAME/texts/memory.txt
du -h $SCRATCH/$NAME/run0/data/*_rhoe0_0.vtk >> $WORK/$NAME/texts/memory.txt
du    $SCRATCH/$NAME/run0/data/*_rhoe0_0.vtk >> $WORK/$NAME/texts/memory.txt
echo "# PARTICLES MEMORY [KB]:" >> $WORK/$NAME/texts/memory.txt
du -h $SCRATCH/$NAME/run0/data/restart0.hdf >> $WORK/$NAME/texts/memory.txt
du    $SCRATCH/$NAME/run0/data/restart0.hdf >> $WORK/$NAME/texts/memory.txt

#-- CREATE 2D CUTS
mkdir $WORK/$NAME/data1
for irun in $(seq 0 $nrun)
do
    echo "$irun"
    python -u test_to2D.py $SCRATCH/$NAME/run$irun/data/ > $WORK/$NAME/texts/RUN_to2D_$irun.out & 
done
wait

#-- MOVE 2D CUTS
for irun in $(seq 0 $nrun)
do
    mv $SCRATCH/$NAME/run$irun/data/*dp*.vtk $SCRATCH/$NAME/run$irun/data/*eq*.vtk $WORK/$NAME/data1
done
cp $SCRATCH/$NAME/run0/data/SimulationData.txt $WORK/$NAME/data1

#-- 3D SELECTION
mkdir $WORK/$NAME/data2
for irun in $(seq 0 $nrun)
do
    for value in {1..9}
    do
        cp $SCRATCH/$NAME/run$irun/data/*_"$value"0.vtk $WORK/$NAME/data2
    done
done
cp $SCRATCH/$NAME/run0/data/*_0.vtk $WORK/$NAME/data2
cp $SCRATCH/$NAME/run0/data/SimulationData.txt $WORK/$NAME/data2

#-- SELECTION PARTICLES
mkdir $WORK/$NAME/dataP
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'Mariner'
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'Bepi'
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'nose'
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'tail'
cat $WORK/$NAME/texts/output_dataP_*.txt > $WORK/$NAME/texts/output_dataP.txt
for irun in $(seq 0 $nrun)
do
    while read line; do 
        cp $SCRATCH/$NAME/run$irun/data/restart"$line".hdf $WORK/$NAME/dataP/restart"$line"_"$irun".hdf 
    done < $WORK/$NAME/texts/output_dataP.txt
done
cp $SCRATCH/$NAME/run0/data/settings.hdf $WORK/$NAME/dataP

#-- PLOT 2D CUTS
python plot_2Dcuts_nBJ.py $WORK/$NAME/data1 $WORK/$NAME/images
python plot_2Dcuts_Pij.py $WORK/$NAME/data1 $WORK/$NAME/images

#-- CREATE 2D CUTS TnV
python -u test_TnV.py $WORK/$NAME/data1 > $WORK/$NAME/texts/RUN_TnV.out &
wait

#-- PLOT 2D CUTS TnV
python plot_2Dcuts_TnV.py $WORK/$NAME/data1 $WORK/$NAME/images

