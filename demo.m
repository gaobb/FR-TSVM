% This is a demo for running linear and nonlinear FRTSVM
%  Author: Bin-Bin Gao 
%  Email: csgaobb@gmail.com
%  July 5, 2017
addpath(genpath('data'));
clc
clear
close all
% load data
load Ripley.mat

%%  Linear FR-TSVM
Parameter.ker='linear';
Parameter.CC=0.25  ;
Parameter.CR= 2^-8;
Parameter.v=10;
Parameter.algorithm='cd';    
Parameter.showplots= 0;
% Training
frtsvm_struct = frtsvmtrain(traindata,trainlabel,Parameter);
fprintf('Num_sv: %d, Time: %.4f s\n', frtsvm_struct.nv, frtsvm_struct.time)
% Testing 
acc = frtsvmclass(frtsvm_struct,testdata,testlabel);


% Nonlinear FR-TSVM
Parameter.ker='rbf';
Parameter.CC=8;
Parameter.CR=1;
Parameter.p1=0.2;
Parameter.v=10;
Parameter.algorithm='cd';   
Parameter.showplots= 1;
% Training
frtsvm_struct = frtsvmtrain(traindata,trainlabel,Parameter);
fprintf('Num_sv: %d, Time: %.4f s\n', frtsvm_struct.nv, frtsvm_struct.time)
% Testing 
acc= frtsvmclass(frtsvm_struct,testdata,testlabel);



