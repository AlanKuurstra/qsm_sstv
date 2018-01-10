% Susceptibility Weighted Imaging and Mapping (SWIM)
%   x = SWIM(threshold, erosion, varargin)
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
%   Haacke et al. JMRI 2010;32(3):663-76
%
%   Created by Tian Liu in 2013
%   Modified by Saifeng Liu on Jan 4, 2014
function x = SWIM(thre, erosion_option, varargin)

[lambda iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir irls smv radius data_weighting gradient_weighting] = parse_QSM_input(varargin{:});

%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
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

g = D;
cone = abs(D)<thre;
g(cone) = sign(g(cone)).*thre;
invg = 1./g;
invg(isnan(invg)) =0;
invg(isinf(invg)) =0;

Z0 = sqrt( (Bx.^2+By.^2)/2);
d = abs(abs(Bz)-Z0);
positive_Za = sqrt((1-3*thre)*(Bx.^2+By.^2)/(2+3*thre));
negative_Za = sqrt((1+3*thre)*(Bx.^2+By.^2)/(2-3*thre));

negative_region = (D<0).*d./abs(negative_Za - Z0);
positive_region = (D>0).*d./abs(positive_Za - Z0);
d = negative_region+positive_region;
d(d>1) = 1;

if erosion_option==1
    %erosion
    [cx cy cz]=ndgrid(-3:3);
    se=(sqrt(cx.^2+cy.^2+cz.^2)<=3);
    Mask=imerode(Mask,se);  
end    
RDF = RDF.*Mask;

invg(isnan(invg)) = 0;    invg(isinf(invg)) = 0;
d(isnan(d)) = 0;    d(isinf(d)) = 0;
X = fftshift(fftn(RDF)).*invg.*d.^2;
x = ifftn(fftshift(X))/(2*pi*delta_TE*CF)*1e6;
