% Total Variation with Split Bregmann (TVSB)
%   x = TVSB(varargin)
%   Adapted from Berkin Bilgic at
%   http://onlinelibrary.wiley.com/store/10.1002/mrm.25029/asset/supinfo/mrm25029-sup-0001-suppinfo.zip?v=1&s=3eec708ad38652e1ba3029ac7a4125db4252f189
%
%   output
%   x - the susceptibility distribution 
%   
%   input
%   default RDF.mat is in current folder.  
%
%   When using the code, please cite 
%   Bilgic et al. MRM 2013 epub
%
%   Adapted from Berkin Bilgic at
%   http://onlinelibrary.wiley.com/store/10.1002/mrm.25029/asset/supinfo/mrm25029-sup-0001-suppinfo.zip?v=1&s=3eec708ad38652e1ba3029ac7a4125db4252f189
function x = TVSB(varargin)
[lambda iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir irls smv radius data_weighting gradient_weighting] = parse_QSM_input(varargin{:});
% lambda = .01;     % L1 penalty

mu = .025;         % Gradient consistency => pick from L2-closed form recon
                    % since the first iteration gives L2 recon
threshold = lambda/mu;

% lambda = threshold*mu;
N = matrix_size;
[k2,k1,k3] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);

fdx = 1 - exp(-2*pi*1i*k1/N(1));
fdy = 1 - exp(-2*pi*1i*k2/N(2));
fdz = 1 - exp(-2*pi*1i*k3/N(3));

cfdx = conj(fdx);           cfdy = conj(fdy);          cfdz = conj(fdz);

E2 = abs(fdx).^2 + abs(fdy).^2 + abs(fdz).^2;
D = dipole_kernel(matrix_size, voxel_size, B0_dir);
D2 = abs(D).^2;

SB_reg = 1 ./ (eps + D2 + mu * E2);
pad_size = [0 0 0];
D = dipole_kernel(matrix_size, voxel_size, B0_dir); %corner
nfm_Sharp_lunwrap = RDF.*Mask;
mask_sharp = Mask;

cfdx = conj(fdx);           cfdy = conj(fdy);          cfdz = conj(fdz);

DFy = conj(D) .* fftn(nfm_Sharp_lunwrap);

SB_reg = 1 ./ (eps + D2 + mu * E2);

vx = zeros(N);          vy = zeros(N);          vz = zeros(N);
nx = zeros(N);          ny = zeros(N);          nz = zeros(N);
Fu = zeros(N);


tic
for t = 1:20
    
    Fu_prev = Fu;
    
    Fu = ( DFy + mu * (cfdx.*fftn(vx - nx) + cfdy.*fftn(vy - ny) + cfdz.*fftn(vz - nz)) ) .* SB_reg;
    
    Rxu = ifftn(fdx .*  Fu);    Ryu = ifftn(fdy .*  Fu);    Rzu = ifftn(fdz .*  Fu);
    
    rox = Rxu + nx;    roy = Ryu + ny;    roz = Rzu + nz;
    
    vx = max(abs(rox) - threshold, 0) .* sign(rox);
    vy = max(abs(roy) - threshold, 0) .* sign(roy);
    vz = max(abs(roz) - threshold, 0) .* sign(roz);
    
    nx = rox - vx;     ny = roy - vy;     nz = roz - vz;
    
    res_change = 100 * norm(Fu(:) - Fu_prev(:)) / norm(Fu(:));
    disp(['Iteration  ', num2str(t), '  ->  Change in Chi: ', num2str(res_change), ' %'])
    
    if res_change < 1
        break
    end
    
end
toc


x = ifftn(Fu) .* mask_sharp/(2.0*pi*delta_TE*CF*1.0e-6);
% chi_SB = real( chi_sb(1+pad_size(1):end-pad_size(1),1+pad_size(2):end-pad_size(2),1+pad_size(3):end-pad_size(3)) ) ;
end