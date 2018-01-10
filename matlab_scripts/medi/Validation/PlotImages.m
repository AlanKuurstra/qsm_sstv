clc; clear all;close all;
addpath('../')
load Truth
load ../HEIDI/qsmHeidi
load ../MEDI/qsmMEDI
load ../TVSB/qsmTVSB
load ../SWIM/qsmiSWIM
load ../TKD/qsmTKD
load ../TSVD/qsmTSVD
load ../R2s/qsmR2s
load ../CSC/qsmCSC


qsmR2s{2} = qsmR2s{2}.*qsmMask{2}*1e-2;
qsmR2s{3} = qsmR2s{3}.*qsmMask{3}*3.75e-3;

brimg(qsmR2s,qsmMask,'R2s');
brimg(qsmTruex,qsmMask,'COSMOS');
brimg(qsmHeidi,qsmMask,'HEIDI');
brimg(qsmTKD,qsmMask,'TKD');
brimg(qsmiSWIM,qsmMask,'iSWIM');
brimg(qsmMEDI,qsmMask,'MEDI');
brimg(qsmTVSB,qsmMask,'TVSB');
brimg(qsmTSVD,qsmMask,'TSVD');
brimg(qsmCSC,qsmMask,'CSC');
