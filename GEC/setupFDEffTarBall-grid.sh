#!/bin/bash
echo "unsetup mrb"
unsetup mrb
# setup mrb v4_04_06
echo "resetup mrb"
setup mrb -o

DIRECTORY=FDEff
# we cannot rely on "whoami" in a grid job. We have no idea what the local username will be.
# Use the GRID_USER environment variable instead (set automatically by jobsub).
USERNAME=${GRID_USER}
echo "Got username"
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
<<<<<<< HEAD
echo "unsetup mrb"
unsetup mrb
# setup mrb v4_04_06
echo "resetup mrb"
setup mrb -o
=======
>>>>>>> 98296839050c17235a06fe4995070adc2451cfcd
setup dunetpc v09_22_02 -q e19:debug # this also sets up the ifdh used in run_FDEffTarBall.sh
echo "dunetpc is set"
export WORKDIR=${_CONDOR_JOB_IWD} # if we use the RCDS the localProducts area will be placed in $CONDOR_DIR_INPUT
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo ~`
fi

# _CONDOR_JOB_IWD is /srv
echo "Check work dir _CONDOR_JOB_IWD: ls -l ${_CONDOR_JOB_IWD}"
ls -l ${_CONDOR_JOB_IWD}

echo "See where are we: pwd"
pwd

# Copy the untarred folder and remove the CVMFS linked read-only version
# because later mrbsetenv needs to modify larsoft-v09_22_02 in build* directory
# which can't be done in the CVMFS read-only version
echo "copy ${INPUT_TAR_DIR_LOCAL}/${DIRECTORY} as ${DIRECTORY}_copy"
cp -r ${INPUT_TAR_DIR_LOCAL}/${DIRECTORY} ${DIRECTORY}_copy

echo "remove CVMFS version: rm -rf ${DIRECTORY}"
rm -rf ${DIRECTORY}

echo "rename copied version: mv ${DIRECTORY}_copy ${DIRECTORY}"
mv ${DIRECTORY}_copy ${DIRECTORY}

echo "Check work dir again _CONDOR_JOB_IWD: ls -l ${_CONDOR_JOB_IWD}"
ls -l ${_CONDOR_JOB_IWD}

echo "Started set up grid"
# INPUT_TAR_DIR_LOCAL is /srv/no_xfer/0/TRANSFERRED_INPUT_FILES/FDEff
# copy FDEff to local work dir _CONDOR_JOB_IWD: /srv
# now is /srv/FDEff/srcs

echo "Check setup-grid exist and source it: ls -l ${_CONDOR_JOB_IWD}/${DIRECTORY}/localProducts*"
ls -l ${_CONDOR_JOB_IWD}/${DIRECTORY}/localProducts*
source ${_CONDOR_JOB_IWD}/${DIRECTORY}/localProducts*/setup-grid

echo "Finished set up grid, now mrbsetenv"
mrbsetenv
echo "Finish mrbsetenv"
