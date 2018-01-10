function Linearplot(true_X,X,qsmMask,name)
%%%%% simu
true_x = real(true_X{1});
x = real(X{1}); 
load('simu_RDF','voxel_size');
% x = fermi_filter(x,size(x), voxel_size);
Mask = qsmMask{1};
filename = ['Simulated Brain' '-' name 'ROI Linear Regression plot'];
ptitle = name;
localLinearplot(true_x,x,Mask,ptitle,filename);
%%%% Gd phantom
true_x = real(true_X{2});
x = real(X{2}); 
load('phantom_RDF','voxel_size');
% x = fermi_filter(x,size(x), voxel_size);
Mask = qsmMask{2};
filename = ['Gd phantom' '-' name 'ROI Linear Regression plot'];
ptitle = name;
localLinearplot(true_x,x,Mask,ptitle,filename);
%%%% in vivo
true_x = real(true_X{3});
x = real(X{3}); %In vivo brain
load('human_RDF','voxel_size');
% x = fermi_filter(x,size(x), voxel_size);

%%%%%% gray matter
% fid = fopen('ROI.img');
% Mask = fread(fid, inf,'ushort');
% fclose(fid);
% Mask = reshape(Mask, size(x));
%%%%%% all mask %%
Mask = qsmMask{3};

filename = ['COSMOS Brain' '-' name 'ROI Linear Regression plot'];
ptitle = name;
localLinearplot(true_x,x,Mask,ptitle,filename);

end

function localLinearplot(true_x,x,Mask,ptitle,name)
ROI_mean = zeros([max(Mask(:)),2] );
ROI_std = zeros([max(Mask(:)),2]);
for r= 1:double(max(Mask(:)));
    ROI_mean(r,1) = mean(true_x(Mask(:)==r));
    ROI_mean(r,2) = mean(x(Mask(:)==r));
    ROI_std(r,1) = std(true_x(Mask(:)==r));
    ROI_std(r,2) = std(x(Mask(:)==r));
end
figure('Visible','off');
% figure;
hold on;

% ROI based
% errorbarxy(ROI_mean(:,1), ROI_mean(:,2), ROI_std(:,1), ROI_std(:,2),{'r','b','b'})
% fitcoef = polyfit(ROI_mean(:,1),ROI_mean(:,2),1);
% R2 = corrcoef(ROI_mean(:,1), ROI_mean(:,2)).^2;R2 = R2(1,2);

% voxel based

plot(true_x(Mask(:)>0),x(Mask(:)>0),'.','markersize',1);
fitcoef = polyfit(true_x(Mask(:)>0),x(Mask(:)>0),1);
R2 = corrcoef(true_x(Mask(:)>0), x(Mask(:)>0)).^2;R2 = R2(1,2);


lo_lim = -0.05;
hi_lim = 0.25;
lo_lim = -0.15;
hi_lim = 0.95;
lfit = fitcoef(1)*[lo_lim:0.01:hi_lim]+fitcoef(2);
plot([lo_lim:0.01:hi_lim],lfit,'-k');
axis([lo_lim hi_lim lo_lim hi_lim]); axis square
if fitcoef(2)>0
    title([ptitle '; y = ' num2str(fitcoef(1),'%02.2f') '*x + ' num2str(fitcoef(2),'%02.2f') '; R^2 = ' num2str(R2,'%02.2f')],'fontsize',16);
else
    title([ptitle '; y = ' num2str(fitcoef(1),'%02.2f') '*x - ' num2str(-fitcoef(2),'%02.2f') '; R^2 = ' num2str(R2,'%02.2f')],'fontsize',16);
end
%%%%% gray matter
% if ~exist('GM_linear','dir');
%     mkdir('GM_linear');
% end
% saveas (gcf,['GM_linear\', name '.tif']);

%%%%% all mask
if ~exist('All_linear','dir');
    mkdir('All_linear');
end
saveas (gcf,['All_linear\', name '.tif']);



end