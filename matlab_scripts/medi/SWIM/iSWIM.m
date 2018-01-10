% iterative Susceptibility Weighted Imaging and Mapping (iSWIM)
%   x = iSWIM(threshold, erosion, varargin)
%
%   output
%   x - the susceptibility distribution 
%   
%   input
%   threshold - dipole truncaiton level
%   erosion - further erode the mask
%   default RDF.mat is in current folder.  
%
%   When using the code, please cite 
%   Tang et al. MRM 2013;69(5):1396-407
%   Haacke et al. JMRI 2010;32(3):663-76
%
%   Created by Tian Liu in 2013
%   Modified by Saifeng Liu on Jan 4, 2014

function x = iSWIM(thre,erosion_option,varargin)

[lambda iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir irls smv radius data_weighting gradient_weighting] = parse_QSM_input(varargin{:});

%%%%%%%%%%%%%%% cone generation %%%%%%%%%%%%%%
Bz_dir = B0_dir;
if (Bz_dir == 1)
    Bz_dir = [1 0 0 ]';
elseif (Bz_dir == 2)
    Bz_dir = [0 1 0 ]';
elseif (Bz_dir==3)
    Bz_dir = [0 0 1]';
end

[Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
    -matrix_size(1)/2:(matrix_size(1)/2-1),...
    -matrix_size(3)/2:(matrix_size(3)/2-1));
    X = X/(matrix_size(1)*voxel_size(1));
    Y = Y/(matrix_size(2)*voxel_size(2));
    Z = Z/(matrix_size(3)*voxel_size(3));
if (reshape(Bz_dir,[1 3])~=[0 0 1])
    By_dir = cross(Bz_dir,[0 0 1]);
else
    By_dir = cross(Bz_dir,[0 1 0]);
end
Bx_dir = cross(By_dir, Bz_dir);

Bx_dir = Bx_dir/norm(Bx_dir);
By_dir = By_dir/norm(By_dir);
Bz_dir = Bz_dir/norm(Bz_dir);
Bx = ( X*Bx_dir(1) + Y*Bx_dir(2) + Z*Bx_dir(3) );
By = ( X*By_dir(1) + Y*By_dir(2) + Z*By_dir(3) );
Bz = ( X*Bz_dir(1) + Y*Bz_dir(2) + Z*Bz_dir(3) );
D = 1/3-  Bz.^2./(Bx.^2+By.^2+Bz.^2);
D(isnan(D)) = 0;

cone = abs(D)<thre;
epsilon = 0.001;
if erosion_option==1
    %erosion of the mask%%
    [cx cy cz]=ndgrid(-3:3);
    se=(sqrt(cx.^2+cy.^2+cz.^2)<=3);
    clear cx cy cz;
    Mask=imerode(Mask,se); 
end
%%%%%% start the iterations %%%%%%%%%%%%%%
x0 = SWIM(thre, erosion_option, varargin{:});
x = x0.*Mask;
M3 = veinMask(real(x).*Mask,0.04,0.1);%%
for i= 1:10
    xvm = M3.*x;
    xvm=edge_preserving_avg(xvm,M3,4,voxel_size(1),voxel_size(2),voxel_size(3));%%
    xvmk = fftshift(fftn(xvm));
    xnk = fftshift(fftn(x)).*(~cone) + xvmk.*cone;
    xn = ifftn(fftshift(xnk));
    if sqrt(sum((xn(Mask==1)-x(Mask==1)).^2)/numel(x(Mask==1)))<epsilon; %%
        disp([num2str(i) ' iterations needed']);
        break
    end
    x = xn;
end
x = real(x).*Mask;%%
