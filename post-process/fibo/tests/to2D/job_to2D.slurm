#!/bin/bash
#MSUB -r 2D # Request name
#MSUB -n 4 # Number of tasks to use
#MSUB -c 6    # Number of threads per task to use
#MSUB -T 10000 # Elapsed time limit in seconds
#MSUB -o RUN_%I.out # Standard output. %I is the job id
#MSUB -e RUN_%I.err # Error output. %I is the job id
#MSUB -q xlarge ##hybrid #skylake #rome #xlarge
#MSUB -m work,scratch
#MSUB -Q normal  ### test  ### queue
#MSUB -A gen12428 ##12622

set -x

ccc_mprun -E'--exclusive' -n 1 python -u test_to2D.py /ccc/work/cont005/gen12428/lavorenf/Mercury_SaeInit/PR0/Normal-bigbox/run0/data &
ccc_mprun -E'--exclusive' -n 1 python -u test_to2D.py /ccc/work/cont005/gen12428/lavorenf/Mercury_SaeInit/PR0/Normal-bigbox/run1/data &

ccc_mprun -E'--exclusive' -n 1 python -u test_to2D.py /ccc/scratch/cont005/gen12428/lavorenf/PR0_minusBz/run0/data &
ccc_mprun -E'--exclusive' -n 1 python -u test_to2D.py /ccc/scratch/cont005/gen12428/lavorenf/PR0_minusBz/run1/data &

wait
