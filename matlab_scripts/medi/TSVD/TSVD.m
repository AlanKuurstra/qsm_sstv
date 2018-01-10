% Truncated Singular Value Decomposition (TSVD)
%   x = TSVD(threshold, varargin)
%
%   output
%   x - the susceptibility distribution 
%   
%   input
%   threshold - dipole truncaiton level
%   default RDF.mat is in current folder.  
%
%   When using the code, please cite 
%   Wharton et al. MRM 2010;63(5):1292-304
%
%   Created by Tian Liu in 2013

function [x] = TSVD(threshold, varargin)

[lambda iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir irls smv radius data_weighting gradient_weighting] = parse_QSM_input(varargin{:});

RDF=RDF.*Mask;

D = dipole_kernel(matrix_size, voxel_size, B0_dir,'kspace');
D(abs(D)<threshold) = inf;
x = real(ifftn(fftn(RDF.*Mask)./(D*2*pi*delta_TE*CF*1e-6)));

