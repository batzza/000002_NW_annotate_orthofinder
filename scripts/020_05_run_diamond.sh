#!/bin/bash
jobid1=$(swarm -f 020_10_make_indexes.swarm --module diamond)
echo $jobid1
jobid2=$(swarm -f 020_15_run_diamond.swarm --module diamond --dependency afterany:$jobid1)
