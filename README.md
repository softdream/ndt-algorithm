## The 2D Normal Distributions Transform( NDT ) Scan Match Method
This repo implemented a 2d NDT algorithm by c++ language **without the pcl libray**. </br>
The only dependencies the code needs are Eigen and Opencv(alternative, only for displaying the scan points).</br>
NDT is a more efficient and less computational scanning matching method compared with ICP. Can be used for 2D laser odometer or SLAM front-end directly.

## How to use
dependencies:
``` shell
Eigen3
Opencv3.4.1 ( Alternative )
``
build:
``` shell
mkdir build
cd build
cmake ..
make
```
use:
in directory: build/bin/
execuate the command:
``` shell
./ndt_test
```
## Test Result

## Reference
[1] Peter Biber, Wolfgang Straaer. The Normal Distributions Transform: A New Approach to Laser Scan 
Matching.
