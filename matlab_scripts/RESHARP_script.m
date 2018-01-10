function [LFS,M1]=RESHARP_script(freqLoc,maskLoc,radius,alpha,LFSFilename)
%addpath /usr/local/medi/BackgroundRemoval/RESHARP
%addpath /usr/local/medi/Common
addpath ./medi/BackgroundRemoval/RESHARP
addpath ./medi/Common
addpath ./NIfTI_20140122
freqObj=load_nii(freqLoc);
maskObj=load_nii(maskLoc);
freq=double(freqObj.img);
mask=double(maskObj.img);
matrix_size=size(freq);
voxel_size=freqObj.hdr.dime.pixdim(2:4);
if isstring(radius) || ischar(radius)
    radius=str2double(radius);
end
if isstring(alpha) || ischar(alpha)
    alpha=str2double(alpha);
end
LFS=RESHARP(freq,mask,matrix_size,voxel_size,radius,alpha);

M1 = int16(SMV(mask, matrix_size, voxel_size, radius)>0.999); %double check that this line matches the mask erosion that happens in RESHARP

tmp=make_nii(LFS);
tmp.hdr=freqObj.hdr;
save_nii(tmp,LFSFilename);

tmp=make_nii(M1);
tmp.hdr=maskObj.hdr;
maskfilename=split(LFSFilename,'.');
maskfilename=maskfilename(1);
maskfilename=char(strcat(maskfilename,'_mask.nii'));
save_nii(tmp,maskfilename);
