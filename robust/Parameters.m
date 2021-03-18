function par = Parameters()
% system
par.L=2;%15;%.0284;
par.g=9.81;
par.In=1/3*par.L^2;
par.Gr=par.g*par.L/2;
% controller
par.tau=.4;%.05;
par.p=20;
par.d=5;
par.a=.8;

par.et=.1;
par.tauh=par.tau*(1+par.et);

par.q=1;

par.h=.0005;
par.r=300;
end