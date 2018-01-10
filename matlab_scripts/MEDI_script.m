p
addpath /usr/local/medi/MEDI
addpath /usr/local/medi/Common
addpath ./NIfTI_20140122

smv=0; %keep background removal separate from MEDI
radius=0;
iFreq=0;

RDFObj=load_nii(RDFLoc);
delta_TE=0.004; %this should be set to something reasonable to ensure optimization error is an appropriate quit condition
RDF=double(RDFObj.img)*delta_TE; %medi wants the frequency to be in radians and multiplied by a uniform echo spacing, make sure your freq map is in rad/s
matrix_size=size(RDF);
voxel_size=RDFObj.hdr.dime.pixdim(2:4);

MaskObj=load_nii(MaskLoc);
Mask=logical(MaskObj.img);

if nargin<6 %if no weights given, then use ones
    N_std=ones(size(RDF))*(.1);
else     
    N_stdObj=load_nii(weightLoc);
    N_std=double(N_stdObj.img);
    N_std=N_std*(sqrt(sum(Mask(:)>0))/norm(N_std(Mask(:)>0)));
end

iMagObj=load_nii(iMagLoc);
iMag=double(iMagObj.img);
B0_dir=[0;0;1]; %this is true for every scan so far...should we make it an input for the future?

%delta_TE=1; %we didn't use their complex frequency fitting, so delta_TE is not wrapped up in our frequency map
merit=1; %iteratively update optimization tuning parameters
CF=298060000; %Hz

lambda=str2double(lambda);
    
%waste of space and time, but it's how they designed their function...might
%be worth rewriting their function to not require a saved .mat file
save([pwd,'/RDF.mat'],'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', 'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'merit','smv');
qsm_MEDI = MEDI_linear('lambda',lambda,'zeropad',[0 0 20],'filename',[pwd,'/RDF.mat']);
qsm_MEDI=qsm_MEDI.*Mask;
tmp=make_nii(qsm_MEDI);
tmp.hdr=RDFObj.hdr;
save_nii(tmp,SuscFilename);

