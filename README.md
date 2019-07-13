# Pseudo-Spectral Time-Domain (PSTD) Method for Simulation of Radar Echoes from Large Heterogeneous Domains

This is a MATLAB-based two-dimensional (2D) PSTD full-wave simulator for solving large-scale (e.g. 10-1000 $\lambda$) low-frequency (e.g. HF) electromagnetic scattering problems with the application of radar sounding of planetary subsurfaces.

Developed by Yang Lei, Darmindra Arumugam, Mark Haynes

Copyright (C) 2019 California Institute of Technology.  Government Sponsorship Acknowledged.

Citation: https://github.com/leiyangleon/PSTD


## Instructions:

* Put the PSTD folder under MATLAB search path
* Run "PSTD.m" from the PSTD folder and refer to it for specific meanings of the individual input and output parameters.
* Refer to and modify "dielectric_scene.m" for accommodating the customized dielectric scene.
* Modify "application_example.m" for individual post-processing applications based on output from running PSTD.m. The current version ingests the far-field output from the PSTD run in the mode of active focusing with user-specified receiver locations, and then reproduces the focused radargram using the far-field output.
* Please also refer to the technical draft associated with this solver (by Lei et al.) for detail of the domain design (e.g. definition of the origin of the coordinates).
* Screenshots of the animation saved as a short video file, named "video.mp4" in the PSTD folder.
