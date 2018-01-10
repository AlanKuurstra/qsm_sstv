function out = MS_L2_QSM(params)

alpha = params.alpha;
N = params.N;
D = params.D;

% For gradient
[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);
E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));
EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;

DtD = abs(D).^2;

tic
out.x = real( ifftn(conj(D) .* fftn(params.phase_unwrap) ./ (DtD + alpha * EtE)));
toc

end