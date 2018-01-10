% Morphology Enabled Dipole Inversion (MEDI) with a linear data fidelity term
%   [x, cost_reg_history, cost_data_history] = MEDI_linear(varargin)
%
%   output
%   x - the susceptibility distribution 
%   cost_reg_history - the cost of the regularization term
%   cost_data_history - the cost of the data fidelity term
%   
%   input
%   default RDF.mat is in current folder.  
%   MEDI_Linear('lambda',lam,...) - lam specifies the regularization parameter
%                               lam is in front of the data fidelity term
%
%   ----optional----   
%   MEDI_Linear('smv', radius,...) - specify the radius for the spherical mean
%                                value operator using differential form
%   MEDI_Linear('merit',...) - turn on model error reduction through iterative
%                          tuning
%   MEDI_Linear('zeropad',padsize,...) - zero pad the matrix by padsize
%
%   When using the code, please cite 
%   J. Liu et al. Neuroimage 2012;59(3):2560-8.
%   T. Liu et al. MRM 2011;66(3):777-83
%   de Rochefort et al. MRM 2010;63(1):194-206
%
%   Adapted from Ildar Khalidov
%   Modified by Tian Liu on 2011.02.01
%   Modified by Tian Liu and Shuai Wang on 2011.03.15
%   Modified by Tian Liu and Shuai Wang on 2011.03.28 add voxel_size in grad and div
%   Modified by Tian Liu on 2013.07.24
%   Last modified by Tian Liu on 2014.12.24

function [x, cost_history, cost_L2_history] = MEDI_linear(varargin)

[lambda iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir merit smv radius data_weighting gradient_weighting Debug_Mode] = parse_QSM_input(varargin{:});

%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
cg_max_iter = 100;
cg_tol = 0.01;
max_iter = 10;
tol_norm_ratio = 0.1;
data_weighting_mode = data_weighting;
gradient_weighting_mode = gradient_weighting;
grad = @cgrad;
div = @cdiv;
iter = 0;

tempn = N_std.*Mask;
if (smv)
    S = SMV_kernel(matrix_size, voxel_size,radius);
    D=S.*dipole_kernel(matrix_size, voxel_size, B0_dir);
    Mask = SMV(Mask, matrix_size,voxel_size, radius)>0.999;
    RDF = iFreq - SMV(iFreq, matrix_size, voxel_size, radius);
    RDF = RDF.*Mask;
    tempn = sqrt(SMV(tempn.^2, matrix_size, voxel_size, radius)+tempn.^2);
else
    D=dipole_kernel(matrix_size, voxel_size, B0_dir);
end

w = dataterm_mask(data_weighting_mode, tempn, Mask);
wG = gradient_mask(gradient_weighting_mode, iMag, Mask, grad, voxel_size,0.9);

res_norm_ratio = Inf;
cost_L2_history = zeros(1,max_iter);
cost_history = zeros(1,max_iter);

e=0.000001; %a very small number to avoid /0
precision_flag = logspace(-1,-16,16);
x = zeros(matrix_size);
while (res_norm_ratio>tol_norm_ratio)&&(iter<max_iter)
tic
    iter=iter+1;
    Vr = 1./sqrt(abs(wG.*grad(real(x),voxel_size)).^2+e);
    A =  @(dx) (1/lambda)*div(wG.*(Vr.*(wG.*grad(real(dx),voxel_size))),voxel_size) + 2*real(ifftn(D.*fftn(w.*w.*real(ifftn(D.*fftn(dx))))));       
    b = A(x) - 2*real(ifftn(D.*fftn(w.*w.*RDF)));

    dx = cgsolve(A, -b, cg_tol, cg_max_iter, 0);
    res_norm_ratio = norm(dx(:))/norm(x(:));
    x = x + dx;

    wres=w.*(real(ifftn(D.*fftn(x))) - RDF);

    cost_L2_history(iter) = norm(wres(:),2);
    cost=abs(wG.*grad(x));
    cost_history(iter) = sum(cost(:));

    
    if merit
        a = wres(Mask==1);
        ma = mean(a); factor = std(a)*5;
        wres = wres-ma;
        wres = abs(wres)/factor;
        wres(wres<1) = 1;
        N_std(Mask==1) = N_std(Mask==1).*wres(Mask==1).^2;
        w(Mask==1) = w(Mask==1)./wres(Mask==1).^2;
    end
    
    fprintf('iter: %d; res_norm_ratio:%8.4f; cost_L2:%8.4f; cost:%8.4f.\n',iter, res_norm_ratio,cost_L2_history(iter), cost_history(iter));
toc
    
end



%convert x to ppm
x = x/(2*pi*delta_TE*CF)*1e6;

if (matrix_size0)
    x = x(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    iMag = iMag(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    RDF = RDF(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    Mask = Mask(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    matrix_size = matrix_size0;
end


end





              