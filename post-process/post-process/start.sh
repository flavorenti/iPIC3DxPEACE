#!/bin/bash
##########################################################
############### START POST-PROCESS #######################
##########################################################
export NAME=test_simu
export WORK=/home/flavorenti/Bureau/Global-simulations-Mercury/
export SCRATCH=/home/flavorenti/Bureau/data_simu_iPIC3D/
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
declare -i irun=0
for run in $SCRATCH/$NAME/run*/ ; do
    cp $run/data/ConservedQuantities.txt $WORK/$NAME/texts/ConservedQuantities$irun.txt
    irun+=1
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

#-- 2D CUTS
mkdir $WORK/$NAME/data1
declare -i irun=0
for run in $SCRATCH/$NAME/run*/ ; do
    echo "$run"
    mpirun -n 1 python -u test_to2D.py $run/data/ > $WORK/$NAME/texts/RUN_to2D_$irun.out &
    irun+=1
done
wait
for run in $SCRATCH/$NAME/run*/ ; do
    mv $run/data/*dp*.vtk $run/data/*eq*.vtk $WORK/$NAME/data1
done
cp $SCRATCH/$NAME/run0/data/SimulationData.txt $WORK/$NAME/data1

#-- 3D SELECTION
mkdir $WORK/$NAME/data2
for run in $SCRATCH/$NAME/run*/ ; do
    for value in {1..9}
    do
        cp $run/data/*_"$value"0.vtk $WORK/$NAME/data2
    done
done
cp $SCRATCH/$NAME/run0/data/*_0.vtk $WORK/$NAME/data2
cp $SCRATCH/$NAME/run0/data/SimulationData.txt $WORK/$NAME/data2

#-- SELECTION PARTICLES
declare -i irun=0
mkdir $WORK/$NAME/dataP
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'Mariner'
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'Bepi'
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'nose'
python print_dataP_list.py $SCRATCH/$NAME/run0/data $WORK/$NAME 'tail'
cat $WORK/$NAME/texts/output_dataP_*.txt > $WORK/$NAME/texts/output_dataP.txt
for run in $SCRATCH/$NAME/run*/ ; do
    while read line; do 
        cp $run/data/restart"$line".hdf $WORK/$NAME/dataP/restart"$line"_"$irun".hdf 
    done < $WORK/$NAME/texts/output_dataP.txt
    irun+=1
done
cp $SCRATCH/$NAME/run0/data/settings.hdf $WORK/$NAME/dataP

#-- PLOT 2D CUTS
python plot_2Dcuts_nBJ.py $WORK/$NAME/data1 $WORK/$NAME/images
python plot_2Dcuts_Pij.py $WORK/$NAME/data1 $WORK/$NAME/images

#-- COMPUTE 2D FIELDS TnV
declare -i irun=0
for run in $SCRATCH/$NAME/run*/ ; do
    mpirun -n 1 python -u test_TnV.py $WORK/$NAME/data1 > $WORK/$NAME/texts/RUN_TnV_$irun.out &
    irun+=1
done
wait

#-- PLOT 2D CUTS TnV
python plot_2Dcuts_TnV.py $WORK/$NAME/data1 $WORK/$NAME/images
