## Processing
The algorithm assumes no phase singularities in the data.  At CFMM, this is accomplished by reconstructing the raw ME-GRE data using an in-house SVD algorithm based on [Klassen](https://cds.ismrm.org/protected/13MProceedings/files/3739.PDF) and [Walsh](https://onlinelibrary.wiley.com/doi/full/10.1002/%28SICI%291522-2594%28200005%2943%3A5%3C682%3A%3AAID-MRM10%3E3.0.CO%3B2-G?sid=nlm%3Apubmed), which gives the least squares best estimate of the magnetization and avoids phase singularities. 

QSM processing was performed as follows. Spatial phase unwrapping was accomplished using a [3D best path algorithm](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-46-26-6623). The frequency at each voxel was then estimated by weighted least squares; each phase echo was weighted by the local SNR in the corresponding T2* image. Finally, background removal and dipole inversion were performed simultaneously using a [single-step QSM algorithm](http://martinos.org/~berkin/Chatnuntawech_2016_NMR_in_Biomed.pdf). 

R2* calculation was estimated using nonlinear least squares. Specifically, python's sipy.optimize.leastsq was used which is a wrapper for [MINPACK's implementation of the Levenberg–Marquardt algorithm](http://netlib.org/minpack/lmder.f). In the signal model ```(A + 1j * B) * np.exp(-TE * R2star + 1j * f * TE)```, ```TE``` is assumed known, while ```R2star```, frequency ```f```, and the compex amplitude ```(A + 1j * B)``` are the unknowns to be estimated. To improve stability, frequency is initialized to the value estimated during qsm processing.

## Run Instructions
[run_instructions.md](run_instructions.md)

## Output
**filename-fieldmap**  
Unwrapped fieldmap in units rad/s without background removal.  
**filename-QSM**  
The susceptibility map in units ppb.  
Note: due to nature of the dipole inversion problem, susceptibility images have an arbitrary offset.  These images have the offset set such that mean value inside the brain is zero.  If you plan to compare values between two different images, it is necessary to segment a structure of known susceptibility and set the offset such that the structure obtains that susceptibility (eg. you can segment the ventricles, take the mean susceptibility value inside the ventricles, and subtract that from your QSM image. This applies an offset that forces the mean ventricle value to be zero.)  
**filename-QSM_brainMask**  
Where the algorithm believes the brain edge is.  
**filename-QSM_noiseMask**  
Voxels that were excluded from the data due to large variance. Usually these are voxels where the phase could not properly be unwrapped, for example near high susceptibility structures like vessels. To reduce the number of voxels masked, increase the value of the input ```scnd_diff_reliability_thresh_noise```.
**filename-R2star**  
The estimated R2star in units Hz from fitting to the complex data.  
**filename-R2star_fit**  
A map of how well different voxels were fit.  
**filename-R2star_negMask**  
Due to noise, R2star values close to zero can sometimes be estimated as a negative value during the fit. Since this is physically impossible, these voxels are set to be 0 but their locations are recorded in this map.  
