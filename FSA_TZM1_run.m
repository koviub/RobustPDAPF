% clear all
% clc
%% Multi-Dimensional Bisection Method
% parameter dimension : 1
% co-dimension (number of equations): 1
addpath('C:\Users\Kovacs Balazs\Desktop\MDBM-Matlab-master\code_folder')
par=Parameters1();
% MDBM
ax=[];
ax(1).val=linspace(0,1,41);  % tau
ax(2).val=linspace(0,10,41);  % L
% ax(3).val=linspace(-1,-.8,5); % eps_tau
bound_function_name='FSA_TZM1';
% if par.et<0
% ax(3).val=linspace(0,50,21);  % om
% bound_function_name='FSA_DZMwC';
% end

par.interpolationorder=3;
Niteration=3;
figure(102)
cmap=colormap;
hold on;grid on;box on;
k=1;
ft=fittype('a1*x.^2');
for eps=linspace(-.9,1,length(cmap))
    par.eb=eps; % 0-->1
    mdbm_sol=mdbm(ax,bound_function_name,Niteration,[],par);
    
Tau=mdbm_sol.posinterp(1,:);
LC=mdbm_sol.posinterp(2,:);    
    
    if eps~=0
    % if wrong fit check exclude data
     fitted=fit(Tau.',LC.',ft,'startpoint',1,'exclude',LC<fitted.a1/4);
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
 plot(x(1:end-1),y(1:end-1),'color',cmap(k,:))
   % set(graphHandle,'EdgeColor',cmap(k,:));
   k=k+1;
end

plot(0:.1:1,3/4*9.81*(0:.1:1).^2,'k')

xlabel('$\tau$[s]','interpreter','latex')
ylabel('$L_{\rm crit}$[m]','interpreter','latex')

figure('position',[10 10 400 400])
hold on;grid on;box on;

plot(valEps,valA1./(3/4*9.81),'r')
xlabel('$\epsilon$','interpreter','latex')
ylabel('$a_1/a_0$','interpreter','latex')
