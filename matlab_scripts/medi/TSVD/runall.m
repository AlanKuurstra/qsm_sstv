clear all;
for i = 1:1;
    tic;runme1; toc;
end

for i = 1:1;
    tic;runme2;toc
end

for i = 1:1;
    tic;runme3;toc;
end
qsm = struct([]);
qsm{1} = qsmsimulation;
qsm{2} = qsmgdphantom;
qsm{3} = qsminvivo;
qsmTSVD = qsm;
save qsmTSVD qsmTSVD