% Truncated K-space Divison (TKD)
%   x = TKD(threshold, varargin)
%
%   output
%   x - the susceptibility distribution 
%   
%   input
%   threshold - dipole truncaiton level
%   default RDF.mat is in current folder.  
%
%   When using the code, please cite 
%   Shmueli et al. MRM 2009;62(6):1510-1522
%
%   Created by Tian Liu in 2013
%   Modified by Karin Shmueli on Jan 3, 2014

function x = TKD(threshold, varargin)

[lambda iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir irls smv radius data_weighting gradient_weighting] = parse_QSM_input(varargin{:});

D = dipole_kernelKS(matrix_size, voxel_size, B0_dir);

threshold = 1/threshold;
D1 = 1/D;
D1((D1>threshold)) = threshold;
D1((1/D<-threshold)) = -threshold;

X = fftn(RDF.*Mask).*D1; 

x = real(ifftn(X)/(2.0*pi*delta_TE*CF*1.0e-6)); 


    