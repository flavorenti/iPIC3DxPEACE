#!/bin/bash
#MSUB -r TnV # Request name
#MSUB -n 2 # Number of tasks to use
#MSUB -c 4    # Number of threads per task to use
#MSUB -T 10000 # Elapsed time limit in seconds
#MSUB -o RUN_%I.out # Standard output. %I is the job id
#MSUB -e RUN_%I.err # Error output. %I is the job id
#MSUB -q xlarge ##hybrid #skylake #rome #xlarge
#MSUB -m work,scratch
#MSUB -Q normal  ### test  ### queue
#MSUB -A gen12428 ##12622

ml python

set -x

ccc_mprun -E'--exclusive' -n 1 python -u test_TnV.py /ccc/work/cont005/gen12428/lavorenf/Mercury_SaeInit/PR0/Normal-bigbox/run0/data1 &
ccc_mprun -E'--exclusive' -n 1 python -u test_TnV.py /ccc/work/cont005/gen12428/lavorenf/Mercury_SaeInit/PR0/Normal-bigbox/run0/data2 &

wait
