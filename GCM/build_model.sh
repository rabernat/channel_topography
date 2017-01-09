#!/bin/bash

MITgcm=/u/rpaberna/MITgcm
module purge
module load comp-intel/11.1.046  mpi/mpt.1.25 netcdf/3.6.0/intel
cd build/


$MITgcm/tools/genmake2 -mods ../code -of $MITgcm/tools/build_options/linux_amd64_ifort+mpi_ice_nas -rootdir $MITgcm && \
make depend && \
make && \
echo "Build successful!"

cd ..
