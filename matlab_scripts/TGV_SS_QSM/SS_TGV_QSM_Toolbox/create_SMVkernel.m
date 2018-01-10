function out = create_SMVkernel_revision(phase_total, chi_mask, min_radius, max_radius, step_size_radius, N, voxel_size)

smv_radii = (max_radius:-step_size_radius:min_radius)*min(voxel_size);
disp(['SMV kernel size = ', num2str(smv_radii)])
num_kernel = length(smv_radii);

[Y,X,Z] = meshgrid(-N(2)/2:(N(2)/2-1),-N(1)/2:(N(1)/2-1),-N(3)/2:(N(3)/2-1));

X = X * voxel_size(1);
Y = Y * voxel_size(2);
Z = Z * voxel_size(3);

unrely_tol = 1e-3;      % tol for unreliable boundary voxels 
SMV_kernel = zeros([N, num_kernel]);
mask_Sharp = zeros([N, num_kernel]);
mask_prev = zeros(N);
phase_Sharp = zeros(N);

count = 1;
for k = 1:num_kernel

    SMV = gen_SMVkernel_voxel_scaled( X, Y, Z, smv_radii(k));
    mask_rely = gen_SMVMask_new( SMV, chi_mask, unrely_tol);   
    
    if sum(mask_rely(:)) == 0
        continue
    end
    
    SMV_kernel(:,:,:,count) = SMV;
    mask_Sharp(:,:,:,count) = (mask_rely-mask_prev);
    phase_Sharp = phase_Sharp +  mask_Sharp(:,:,:,count).* ifftn(SMV_kernel(:,:,:,count) .* fftn(phase_total));
    mask_prev = mask_rely;
        
    if count == 1
        Del_Sharp_Inv = SMV;
    end
    count = count + 1;
end
    

out.SMV_kernel = SMV_kernel;
out.SMV_inv_kernel = Del_Sharp_Inv; % For inverse Sharp
out.SMV_mask = mask_Sharp;
out.mask_eval = phase_Sharp~=0;
out.SMV_phase = phase_Sharp;

end


function SMV = gen_SMVkernel_voxel_scaled( X, Y, Z, smv_rad)
  
smv = (X.^2 + Y.^2 + Z.^2) <= smv_rad^2;
smv = smv / sum(smv(:));

smv_kernel = zeros(size(X));
smv_kernel(1+end/2,1+end/2,1+end/2) = 1;
smv_kernel = smv_kernel - smv;

SMV = fftn(fftshift(smv_kernel));

end

function mask_rely = gen_SMVMask_new( SMV, chi_mask, unrely_tol)
  
mask_unrely = ifftn(SMV .* fftn(chi_mask));
mask_unrely = abs(mask_unrely) > unrely_tol;    % mask of unreliable phase estimates

mask_rely = chi_mask .* (chi_mask - mask_unrely);       % mask of reliable phase estimates

end