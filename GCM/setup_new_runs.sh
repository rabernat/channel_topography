#!/bin/bash

refdir="run_taux2000_rb0110_bump"
mydirs="run_taux0125_rb0110_flat run_taux0250_rb0110_bump run_taux0250_rb0110_flat run_taux0500_rb0110_bump run_taux0500_rb0110_flat run_taux1000_rb0110_bump run_taux1000_rb0110_flat run_taux2000_rb0110_bumplong run_taux4000_rb0110_bump run_taux4000_rb0110_flat run_taux8000_rb0110_bump run_taux8000_rb0110_flat"
for d in $mydirs
do
    cd $d && \
    rm -v data.diagnostics && \
    cp -v "../$refdir/PBS_run_script_short" "../$refdir/data.diagnostics" "../$refdir/data.pkg" . && \
    setup_pickup.sh && \
    qsub PBS_run_script_short && \
    cd ..
done


