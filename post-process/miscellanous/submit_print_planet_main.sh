#!/bin/bash

## PR0 --> 0-56
## PR1 --> 0-86
## PR2 --> 0-87

for value in {0..56}
do
    echo $value
    ccc_msub job_print_planet.irene
    sleep 10
    perl -i -pe 's/inow=\K\d+/'$[ $value ]'/ge' print_planet_main.py
    perl -i -pe 's/RUN_PR-\K\d+/'$[ $value ]'/ge' job_print_planet.irene
done

