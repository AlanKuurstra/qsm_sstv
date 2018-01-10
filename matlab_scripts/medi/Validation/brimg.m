
function brimg(X,M, name2)
x = real(X{1}); %Simulated brain
name = 'Simulated Brain';
load('simu_RDF','voxel_size');
Mask = M{1};
slice = [159 156 49];
texthere1 = 'xz';
% x = fermi_filter(x,size(x), voxel_size);
localimg(x, Mask, voxel_size, slice, [-.2 .2], name, name2,[0 90 90],texthere1);

x = real(X{2}); %Gadolinium phantom
name = 'Gd Phantom';
load('phantom_RDF','voxel_size');
Mask = M{2};
slice = [66 65 75];
texthere2 = 'xy';
% x = fermi_filter(x,size(x), voxel_size);
localimg(x, Mask, voxel_size, slice, [-.45 .45], name, name2,[-90 90 90],texthere2);

x = real(X{3}); %In vivo brain
name = 'COSMOS Brain';
load('human_RDF','voxel_size');
Mask = M{3};
slice = [150 181 66]; 
texthere3 = 'xy';
% x = fermi_filter(x,size(x), voxel_size);
localimg(x, Mask, voxel_size, slice, [-.25 .25], name, name2, [-90 90 90],texthere3);

texthere3 = 'yz';
localimg(x, Mask, voxel_size, slice, [-.20 .20], name, name2, [-90 90 90],texthere3);

end

function localimg(x, Mask, voxel_size, slice, wl, name, name2, rotate,texthere)

matrix_size = size(x);
imgyz = squeeze(x(slice(1),:,:).*Mask(slice(1),:,:));
[x1, x2, y1, y2] = localfindc(squeeze(Mask(slice(1),:,:)));
imgyz = imgyz(x1-2:x2+2,y1-2:y2+2,:);


imgxz = squeeze(x(:,slice(2),:).*Mask(:,slice(2),:));
[x1, x2, y1, y2] = localfindc(squeeze(Mask(:,slice(2),:)));
imgxz = imgxz(x1-2:x2+2,y1-2:y2+2,:);


imgxy = squeeze(x(:,:,slice(3)).*Mask(:,:,slice(3)));
[x1, x2, y1, y2] = localfindc(squeeze(Mask(:,:,slice(3))));
imgxy = imgxy(x1-2:x2+2,y1-2:y2+2,:);

matrix_size = [size(imgxy,1) size(imgxy,2) size(imgxz,2)];
    
imgyz = imresize(imgyz,voxel_size([2 3]).*matrix_size([2 3]));
imgxz = imresize(imgxz,voxel_size([1 3]).*matrix_size([1 3]));
imgxy = imresize(imgxy,voxel_size([1 2]).*matrix_size([1 2]));

matrix_size = [size(imgxy) size(imgxz,2)];
dif = mod(matrix_size,[2 2 2]);
matrix_size = matrix_size - dif;
if dif(1)==1
    imgxy = imgxy(1:end-1,:);
    imgxz = imgxz(1:end-1,:);
end
if dif(2)==1
    imgxy = imgxy(:,1:end-1);
    imgyz = imgyz(1:end-1,:);
end
if dif(3)==1
    imgxz = imgxz(:,1:end-1);
    imgyz = imgyz(:,1:end-1);
end
numpad = ([max(matrix_size)-matrix_size(1) max(matrix_size)-matrix_size(2) max(matrix_size)-matrix_size(3)]/2);
imgxy = padarray(imgxy,numpad([1 2]));
imgxz = padarray(imgxz,numpad([1 3]));
imgyz = padarray(imgyz,numpad([2 3]));

imgxy = imrotate(imgxy,rotate(1));
imgxz = imrotate(imgxz,rotate(2));
imgyz = imrotate(imgyz,rotate(3));

img = [imgxy imgxz imgyz];
if strcmp(texthere, 'xy')
    img = [imgxy];
elseif strcmp(texthere, 'yz')
    img = [imgyz];
elseif strcmp(texthere, 'xz')
    img = [imgxz];
end
img = (img-wl(1))/(wl(2)-wl(1));
I = imresize(img,size(img)*2);
I2 = real(I);
if ~exist('QSM_images','dir')
    mkdir('QSM_images');
end
imwrite(double(I2),['QSM_images\' name '-' name2 '-' texthere '.bmp'],'bmp');

end

function [x1, x2, y1, y2] = localfindc(img)
img(img~=0)=1;
I = max(max(img,[],2),[],3);
x1 = max(find(I,1,'first'));
x2 = min(find(I,1,'last'));
I = max(max(img,[],1),[],3);
y1 = max(find(I,1,'first'));
y2 = min(find(I,1,'last'));
end

