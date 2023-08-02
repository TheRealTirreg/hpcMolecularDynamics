#! /bin/sh

echo "Building 08"
echo "========================================================="
cd "../../cmake-build-release"
pwd
ls -la
echo "========================================================="
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON ..
ninja
echo "========================================================="
pwd
ls -la
echo "========================================================="
mpirun -n 4 "milestones/08/08"