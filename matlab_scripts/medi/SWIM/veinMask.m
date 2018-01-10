function M3 = veinMask(x, thre1, thre2)
%Modified by Saifeng Liu on Jan 4, 2014
if nargin<3
    thre2 = 0.25;
end
if nargin<2
    thre1 = 0.07;
end

for z = 1:size(x,3)
    x(:,:,z) = medfilt2(x(:,:,z),[3 3]);
end

x1 = zeros(size(x)); x1(:,:,1:end-2) = x(:,:,3:end);
x2 = zeros(size(x)); x2(:,:,1:end-1) = x(:,:,2:end);
x3 = zeros(size(x)); x3(:,:,1:end) = x(:,:,1:end);
x4 = zeros(size(x)); x4(:,:,2:end) = x(:,:,1:end-1);
x5 = zeros(size(x)); x5(:,:,3:end) = x(:,:,1:end-2);
MIP = max(cat(4,x1,x2,x3,x4,x5),[],4);

M0 = x>thre1;
MP = MIP>thre2;
[x y z]=ndgrid(-3:3);
se=(sqrt(x.^2+y.^2+z.^2)<=3);
M1=bwareaopen(M0,125);%remove small islands
M2=imclose(M1,se);%remove small holes
M2=~bwareaopen(~M2,64);
M3 = M2.*MP;


