#!/bin/bash

MITgcm=/u/rpaberna/MITgcm
module purge
#module load comp-intel/11.1.046  mpi/mpt.1.25 netcdf/3.6.0/intel
module load comp-intel/2011.7.256 mpi-sgi/mpt.2.06a67 netcdf/4.0
cd build_long/

$MITgcm/tools/genmake2 -mods ../code_long -of ../code_long/linux_amd64_ifort+mpi_ice_nas_SB -rootdir $MITgcm && \
make depend && \
make && \
echo "Build successful!"

cd ..
