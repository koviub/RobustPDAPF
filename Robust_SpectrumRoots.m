% clc
clear all

uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY/'ALL''
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();
tau0=par.tau;
p0=par.p;
d0=par.d;
a0=par.a;
L0=par.L;

%% invetistigated domain
x0=0+1i*0;  % origin on the complex plane
Rad=10;
%% modell

[W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,'STRUCTURED');
%% Exact roots

root=TransRoot(Dc,4,x0,Rad);
% figure()
hold on
plot(root,'ro','markersize',6,'markerfacecolor','r')
plot([0 0],[-Rad Rad],'b','linewidth',1.5)
plot([-Rad Rad],[0 0],'k--','linewidth',1)
xlim([-Rad Rad])
ylim([-Rad Rad])

for i=1:50
    
    par=PerturbParams(par,uncertainty,perturbation,[L0,p0,d0,a0,tau0]);
    
    % measure of perturbation
    switch lower(uncertainty)
        case 'delay'
            epsi=abs((par.tau-tau0)/tau0);
        case 'system'
            epsi=abs((par.L-L0)/L0);
        case 'control'
            epsi=norm(abs([(par.p-p0)/p0,(par.d-d0)/d0,(par.a-a0)/a0]));%abs((par.a-a0)/a0);
        case 'all'    
            epsi=norm(abs([(par.L-L0)/L0,(par.p-p0)/p0,(par.d-d0)/d0,(par.a-a0)/a0]));%abs((par.a-a0)/a0);
       
    end
    
    [W,Mx,Dc,w,G]=ModelDefinition(par,uncertainty,'STRUCTURED');
    %% Exact roots
    
    root=TransRoot(Dc,4,x0,Rad);
    
    hold on
    c = epsi*ones(size(root));
    scatter(real(root),imag(root),10,c,'filled','MarkerEdgeColor','none')
    drawnow();
    
end

function par=PerturbParams(par,str0,str1,par0)
L0=par0(1);
p0=par0(2);
d0=par0(3);
a0=par0(4);
tau0=par0(5);
str1=lower(str1);

rszog=rand;
switch lower(str0)
    case 'system'
        switch str1
            case 'complex'
                par.L=L0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*L0;
                
            case 'real'
                par.L=L0+(2*rand-1)*L0;
        end
        
    case 'control'
        switch str1
            case 'complex'
                par.p=p0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*p0;
                par.d=d0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*d0;
                par.a=a0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*a0;
                
            case 'real'
                par.p=p0+(2*rand-1)*p0;
                par.d=d0+(2*rand-1)*d0;
                par.a=a0+(2*rand-1)*a0;
                
        end
        
    case 'delay'
        switch str1
            case 'complex'
                par.tau=tau0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*tau0;
                
            case 'real'
                par.tau=tau0+(2*rand-1)*tau0;
        end
        
    case 'all'
        switch str1
            case 'complex'
                par.p=p0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*p0;
                par.d=d0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*d0;
                par.a=a0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*a0;
                par.L=L0+rand*(cos(2*pi*rszog)+1i*sin(2*pi*rszog))*L0;
            case 'real'
                par.p=p0+(2*rand-1)*p0;
                par.d=d0+(2*rand-1)*d0;
                par.a=a0+7*(2*rand-1)*a0;
                par.L=L0+(2*rand-1)*L0;
        end
end
end