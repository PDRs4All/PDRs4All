#!/usr/bin/env bash
set -e

# -- multiprocessing options --
# _____________________________

# J is the number of processes for stage 1 and 2. The recommended limit, is to make sure you
# have about 10 GB of RAM per process. For the science cluster at ST, with 512 GB RAM, I use
# J=48.
J=1

# Use these if there's too much multithreading. On machines with high core counts, numpy etc can
# sometimes launch a large number of threads. This doesn't give much speedup if multiprocessing
# is used already.
T=8
export MKL_NUM_THREADS=$T
export NUMEXPR_NUM_THREADS=$T
export OMP_NUM_THREADS=$T
export OPENBLAS_NUM_THREADS=$T

# -- environment --
# _________________

# Set CRDS context here. If N is not a number, no context will be set, resulting in the latest
# pmap.
export CRDS_PATH=$HOME/crds_cache
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
# N=1147
N=latest
if [[ $N =~ ^[0-9]{4}$ ]]
then export CRDS_CONTEXT=jwst_$N.pmap
fi

# If you use conda, activate enable its use here and activate the right environment
# eval "$(conda shell.bash hook)"
# conda activate jwstpip
python -c "import sys; print(sys.executable)"

# -- set up directories --
# ________________________

# Specify input directories as recommended in the readme. The script needs absolute paths for
# everything. This is a quirk of the association generator.
HERE=$(pwd)
IN_SCI=$HERE/science
IN_SCII=$HERE/science_imprint
IN_BKG=$HERE/background
IN_BKGI=$HERE/background_imprint

# Modify the output directory here. Default has the pmap number in it (if set).
OUT_PFX=${N}pmap
OUT_SCI=$HERE/$OUT_PFX/science
OUT_SCII=$HERE/$OUT_PFX/science_imprint
OUT_BKG=$HERE/$OUT_PFX/background
OUT_BKGI=$HERE/$OUT_PFX/background_imprint

# -- set up logging --
# --------------------
OUT_LOG=$HERE/$OUT_PFX/log
mkdir -p $OUT_LOG
python -c 'import jwst; print(jwst.__version__)' > $OUT_LOG/versions.txt
python -c 'import crds; print(crds.get_context_name("jwst"))' >> $OUT_LOG/versions.txt

parallel_shorthand () {
    echo $2
    cp jobs_$2.sh $OUT_LOG/jobs_$2.sh # save the jobs too. Good to know the exact commands
    parallel --retries 3 --retry-failed --joblog $OUT_LOG/$2.joblog --progress -j $1 {} ">>"$OUT_LOG/$2_cpu{%}.log '2>&1' :::: jobs_$2.sh
}

# -- run the pipeline --
# ______________________

# background imprint
pipeline_jobs -s 1 -o $OUT_BKGI $IN_BKGI
mv strun_calwebb_detector1_jobs.sh jobs_bkgi_1.sh
parallel_shorthand 1 bkgi_1

# background
pipeline_jobs -s 1 -o $OUT_BKG $IN_BKG
mv strun_calwebb_detector1_jobs.sh jobs_bkg_1.sh
parallel_shorthand $J bkg_1

# science imprint
pipeline_jobs -s 1 -o $OUT_SCII $IN_SCII
mv strun_calwebb_detector1_jobs.sh jobs_scii_1.sh
parallel_shorthand $J scii_1

# science
pipeline_jobs -s 1 -o $OUT_SCI $IN_SCI
mv strun_calwebb_detector1_jobs.sh jobs_sci_1.sh
parallel_shorthand $J sci_1

# NSClean is now run in the pipeline. Results are slightly different for the following reasons:
# - There are more settings
# - The mask is created on-the-fly
# - The subtraction is done at every frame before ramp fit step
# Therefore we keep nirspec_nsclean.bash around for now.

# background stage 2
pipeline_jobs -s 2 -i $OUT_BKGI -o $OUT_BKG $IN_BKG
mv strun_calwebb_spec2_jobs.sh jobs_bkg_2.sh
parallel_shorthand 1 bkg_2

# science stage 2
pipeline_jobs -s 2 -i $OUT_SCII -o $OUT_SCI $IN_SCI
mv strun_calwebb_spec2_jobs.sh jobs_sci_2.sh
parallel_shorthand 4 sci_2

# science stage 3
pipeline_jobs -s 3 --mosaic -b $OUT_BKG -o $OUT_SCI $IN_SCI
mv strun_calwebb_spec3_jobs.sh jobs_sci_3.sh
parallel_shorthand 1 sci_3
