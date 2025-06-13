## Installation
To install the package, the right development tools must be in place. Specifically,
Windows users need [RTools](https://cran.r-project.org/bin/windows/Rtools/), while
Mac users need [Xcode and a Fortran compiler](https://cran.r-project.org/bin/macosx/tools/), 
and type the following commands in a terminal to download and extract an OpenMP runtime:
```
curl -O https://mac.r-project.org/openmp/openmp-17.0.6-darwin20-Release.tar.gz
sudo tar fvxz openmp-17.0.6-darwin20-Release.tar.gz -C /
```
Linux users need a GCC compiler. Most Linux systems in a server environment already have 
a GCC compiler with OpenMP support available. If that is not the case, you can try
```
sudo apt install build-essential
```
To install the pacakge, try:
```R
library(devtools)
install_github("qhengncsu/FastSQR")
```
