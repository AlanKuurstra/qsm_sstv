clc;clear all;

for i = 1:1;
    runme1; 
end

for i = 1:1;
    runme2;
end

for i = 1:1;
    runme3;
end

qsm = struct([]);
qsm{1} = qsmsimulation;
qsm{2} = qsmgdphantom;
qsm{3} = qsminvivo;
qsmMEDI = qsm;
save qsmMEDI qsmMEDI