% clear all
% clc
%% Multi-Dimensional Bisection Method
% parameter dimension : 1
% co-dimension (number of equations): 1
addpath('C:\Users\Kovacs Balazs\Desktop\MDBM-Matlab-master\code_folder')
par=Parameters();
% MDBM
ax=[];
ax(1).val=linspace(0,1,41);  % tau
ax(2).val=linspace(0,10,41);  % L
% ax(3).val=linspace(-1,-.8,5); % eps_tau
bound_function_name='FSA_TZM';
% if par.et<0
% ax(3).val=linspace(0,50,21);  % om
% bound_function_name='FSA_DZMwC';
% end

par.interpolationorder=3;
Niteration=3;
figure(102)
cmap=colormap;
NN=length(cmap);
% eps=[linspace(-1,0,NN/2) linspace(1/32,1,N/2)];

hold on;grid on;box on;
k=1;
ft=fittype('a1*x.^2');
for eps=[linspace(-1,0,NN) linspace(1/32,1,NN)]
    par.et=eps; % 0-->1
    mdbm_sol=mdbm(ax,bound_function_name,Niteration,[],par);
Tau=mdbm_sol.posinterp(1,:);
LC=mdbm_sol.posinterp(2,:);    
    if eps~=0
    % if wrong fit check exclude data
    if eps<-.5
    fitted=fit(Tau.',LC.',ft,'startpoint',1,'exclude',LC<4);
    else
    fitted=fit(Tau.',LC.',ft,'startpoint',1);
    end
    valA1(k)=fitted.a1;
    valEps(k)=eps;
    %         hold on
    %         plot(fitted,valT(:,h),valL(:,h),(valL(:,h)>=exc))
    %         legend('off')
    else
        valA1(k)=0;
        valEps(k)=eps;
    end
    x=sort(Tau).';
    y=(fitted.a1*x.^2);
%     graphHandle=plot_mdbm(mdbm_sol,[1 0 0]);view([0,0,1]);
    plot(x(1:end-1),y(1:end-1),'color',cmap(ceil(k/2),:))
%     set(graphHandle,'EdgeColor',cmap(k,:));
    k=k+1;
end

plot(0:.1:1,3/4*9.81*(0:.1:1).^2,'k')

xlabel('$\tau$[s]','interpreter','latex')
ylabel('$L_{\rm crit}$[m]','interpreter','latex')


[tt,ep]=meshgrid(0:.05:1,valEps);
ep(ep==Inf)=1;
LL=[];
for i=1:size(tt,2)
    for j=1:size(ep,1)
        LL(j,i)=valA1(j)*tt(1,i)^2;
    end
end
figure();
surface(tt,ep,LL,'CData',abs(ep))
zlim([0 14])
ylim([-1 1])
box on
grid on

figure('position',[10 10 400 400])
hold on;grid on;box on;

plot(valEps(1:end-1),valA1(1:end-1)./(3/4*9.81),'r')
set(gca,'xTickLabel',{'0','0.5\tau','\tau','1.5\tau','2\tau'})
xlabel('$\hat\tau$','interpreter','latex')
ylabel('$c(\Delta)/c_0$','interpreter','latex')
