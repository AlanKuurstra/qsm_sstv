clc;clear all;close all;
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


Linearplot(qsmTruex,qsmHeidi,qsmMask,'HEIDI');
Linearplot(qsmTruex,qsmTKD,qsmMask,'TKD');
Linearplot(qsmTruex,qsmiSWIM,qsmMask,'iSWIM');
Linearplot(qsmTruex,qsmMEDI,qsmMask,'MEDI');
Linearplot(qsmTruex,qsmTVSB,qsmMask,'TVSB');
Linearplot(qsmTruex,qsmTSVD,qsmMask,'TSVD');
Linearplot(qsmTruex,qsmCSC,qsmMask,'CSC');
