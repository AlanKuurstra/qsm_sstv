%Edge preserving averaging
%by Saifeng Liu on Jan 4, 2014
function ima_lp=edge_preserving_avg(ima,veinMask,radius,dx,dy,dz)
%dx x dy x dz is the voxel size
dim=size(ima);
[cy,cx,cz]=meshgrid(-radius:radius,-radius:radius,-radius:radius);
rho_sp=zeros(2*radius+1,2*radius+1,2*radius+1);
index_sp=(cy*dy).^2+(cx*dx).^2+(cz*dz).^2<=radius^2;%a small sphere
num=sum(index_sp(:));
rho_sp(index_sp)=1;
rho=zeros(dim);
rho(dim(1)/2+1-radius:dim(1)/2+1+radius,dim(2)/2+1-radius:dim(2)/2+1+radius,dim(3)/2+1-radius:dim(3)/2+1+radius)=rho_sp;
rhofft=fftn(rho);
clear rho;
ima_lp=ifftshift(ifftn(fftn(ima).*rhofft))/num;    
lpMask=ifftshift(ifftn(fftn(veinMask).*rhofft))/num;    
lpMask(veinMask==0)=1;
ima_lp=ima_lp./lpMask;
ima_lp=ima_lp.*veinMask;
end