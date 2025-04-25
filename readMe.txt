This is a package of FID/SE hybrid MRSI reconstruction code, including:
(1) Main reconstruction function: Jresolved_MRSI_FIDSE_jointRecon.m
(2) Support functions in folder "support"
(3) Demo script: demo_Jresolved_MRSI_FIDSE_jointRecon.m

Demo data can be downloaded from http://mri.beckman.illinois.edu/software.html, including:
 - demo_jointRecon.mat: demo data for FID/SE hybrid MRSI reconstruction
 - demo_original_data.mat: demo original acquired data (FID data have been coil-combined to reduce file size)
 - demo_molecular_maps.mat: demo final processed molecular maps

System requirements:
The reconstruction code requires a standard computer with enough RAM to process the FID/SE data.
This package is supported for the Linux operating system. The package has been tested on Linux Ubuntu 16.04 and Linux Ubuntu 18.04.

Software requirements:
The package has been tested on MATLAB R2014a and MATLAB R2023a.

Setup:
All dependencies have been included in the package.
To install, simply extract FIDSE_jointRecon.zip.
The installation should take less than 1 min.

Demo:
To run the demo script, type demo_Jresolved_MRSI_FIDSE_jointRecon in the MATLAB command line.
The results are reconstructed FID/SE MRSI data.
It takes about 45 min for the demo code to finish the MRSI reconstruction on an average workstation.