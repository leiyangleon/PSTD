# PSTD

Pseudo-Spectral Time-Domain Method for Simulation of Radar Echoes from Large Heterogeneous Domains

MATLAB-based 2D PSTD Solver for Radar Sounding Applications

Copyright (C) 2019 Yang Lei, Darmindra Arumugam, Mark Haynes
California Institute of Technology

## Instructions:

* 1. Put the PSTD folder under MATLAB search path
* 2. Run "PSTD.m" from the PSTD folder and refer to it for specific meanings of the individual input and output parameters.
* 3. Refer to and modify "dielectric_scene.m" for accommodating the customized dielectric scene.
* 4. Modify "application_example.m" for individual post-processing applications based on output from running PSTD.m. The current version ingests the far-field output from the PSTD run in the mode of active focusing with user-specified receiver locations, and then reproduces the focused radargram using the far-field output.
* 5. Please also refer to the technical draft associated with this solver (by Lei et al.) for detail of the domain design (e.g. definition of the origin of the coordinates).
* 6. Screenshots of the animation saved as a short video file, named "video.mp4" in the PSTD folder.
