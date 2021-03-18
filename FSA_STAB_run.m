% clear all
% clc
%% Multi-Dimensional Bisection Method
% parameter dimension : 1
% co-dimension (number of equations): 1
addpath('C:\Users\Kovacs Balazs\Desktop\MDBM-Matlab-master\code_folder')


par.L=2;


ax=[];
ax(1).val=linspace(0,200,21);  % p
ax(2).val=linspace(-10,40,21);  % d
ax(3).val=linspace(0,100,21);  % om

figure(104)
hold on
par.interpolationorder=2;
Niteration=5;
bound_function_name='FSA_STAB';
count=1;
for i=0:1/7:1
    subplot(4,2,count);
    par.tau=.25;
    par.et=-1*i; % 0-->1
    mdbm_sol=mdbm(ax,bound_function_name,Niteration,[],par);
    plot_mdbm(mdbm_sol);
    
    view([0 0 1]);
    drawnow();
    % bound_function_name='FSA_STAB2';
    % mdbm_sol2=mdbm(ax,bound_function_name,Niteration,[],par);
    count=count+1;
end






% plot_mdbm(mdbm_sol2);