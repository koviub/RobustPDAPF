% clear all
% clc
%% Multi-Dimensional Bisection Method
% parameter dimension : 1
% co-dimension (number of equations): 1
addpath('C:\Users\Kovacs Balazs\Desktop\MDBM-Matlab-master\code_folder')
par=Parameters();
par.et=-.5; % 0-->1
% MDBM
ax=[];
ax(1).val=linspace(0,1,21);  % tau
ax(2).val=linspace(0,10,21);  % L
ax(3).val=linspace(10,50,31);  % om
bound_function_name='FSA_DZMwC';


par.interpolationorder=2;
Niteration=3;
mdbm_sol=mdbm(ax,bound_function_name,Niteration,[],par);


figure(103)
plot(0:0.1:1,3*9.81/4*(0:.1:1).^2,'r')
hold on
plot_mdbm(mdbm_sol);view([0,0,1])