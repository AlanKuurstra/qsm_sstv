function S=reg_term(xx,wG,mask,dimension)


if nargin<2
    w2=1;
end
if nargin<3
    dimension=[1 1 1];
end
       
x = zeros(size(mask));
x(mask(:) == 1) = xx(1:end);
x(mask(:) == 0) = 0;

S = cdiv(wG.*wG.*cgrad(real(x),dimension), dimension);
S = S(mask(:) == 1);


