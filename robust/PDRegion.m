function [ps,pe,ds,de]=PDRegion(par,num)

if nargin<2
    dr=false;
else
    dr=true;
end

addpath('C:\Users\Kovacs Balazs\Desktop\MDBM-Matlab-master\code_folder')

b=3*9.81/2/par.L;
tau=par.tau;
et=par.et;
om=0.01:.01:100;


% FSA
p=b.^(1/2).*(b+om.^2).*(b.^(1/2).*om.*(3+(-2).*cos(et.*om.* ...
  tau)+4.*cosh(b.^(1/2).*(1+et).*tau).*sin((1/2).*et.*om.*tau) ...
  .*sin((1/2).*(2+et).*om.*tau))+2.*(b+(-1).*om.^2).*cos((1/2) ...
  .*(2+et).*om.*tau).*sin((1/2).*et.*om.*tau).*sinh(b.^(1/2).* ...
  (1+et).*tau)).^(-1).*(om+2.*sin((1/2).*et.*om.*tau).*(om.* ...
  cosh(b.^(1/2).*(1+et).*tau).*sin((1/2).*(2+et).*om.*tau)+ ...
  b.^(1/2).*cos((1/2).*(2+et).*om.*tau).*sinh(b.^(1/2).*(1+et) ...
  .*tau)));
d=(-2).*(b+om.^2).*sin((1/2).*et.*om.*tau).*(b.^(1/2).*om.*(3+ ...
  (-2).*cos(et.*om.*tau)+4.*cosh(b.^(1/2).*(1+et).*tau).*sin(( ...
  1/2).*et.*om.*tau).*sin((1/2).*(2+et).*om.*tau))+2.*(b+(-1) ...
  .*om.^2).*cos((1/2).*(2+et).*om.*tau).*sin((1/2).*et.*om.* ...
  tau).*sinh(b.^(1/2).*(1+et).*tau)).^(-1).*(b.^(1/2).*cos(( ...
  1/2).*(2+et).*om.*tau).*cosh(b.^(1/2).*(1+et).*tau)+om.*sin( ...
  (1/2).*(2+et).*om.*tau).*sinh(b.^(1/2).*(1+et).*tau));

%% Multi-Dimensional Bisection Method
% parameter dimension : 1
% co-dimension (number of equations): 1

ax=[];
ax(1).val=linspace(0,100,51);  % om
par.interpolationorder=2;
Niteration=5;
bound_function_name='PDmax_fval';
mdbm_sol=mdbm(ax,bound_function_name,Niteration,[],par);
frek_end=mdbm_sol.posinterp;
index=find(om>frek_end(1));

de1=d(index(1));
[pk,loc]=findpeaks(d);
de0=pk(1);

if index(1)<loc(1)
    de=de1;
else
    de=de0;
end
if de>4*b
de=4*b;
end

ps=p(1);
[pk,~]=findpeaks(p);
pe=pk(1);
if pe>4*b
pe=4*b;
end

[dk,~]=findpeaks(-d);
ds=-dk(1);
if dr
% figure(num);
hold on;grid on
plot(p,d,'b')
plot([b b],[-de de]*2,'r')
plot([-1 2]*pe,[0 0],'k:')
plot([0 0],[-1,2]*de,'k:')
rectangle('position',[ps ds abs(pe-ps) abs(de-ds)])
xlim([ps/3-1 pe+1])
ylim([ds-1 de+1])
end

end
