clear all
% clc

uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY/'ALL''
structure    = 'STRUCTURED'; % 'STRUCTURED'/'UNSTRUCTURED'
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();

%% modell

[W,Mx,Dc,w,w1,w2,G]=ModelDefinition1(par,uncertainty,structure);

%% invetistigated domain
x0=0+1i*0;                          % origin on the complex plane
Rad=10;                         % radius around x0

re=-Rad:Rad/200:Rad/3;
im=-Rad:Rad/200:Rad;
[R,I]=meshgrid(re,im);
lambda=R+1i*I;
%% Robustness

r=SpectralRadius(W,lambda,perturbation);

figure();
% surface(R,I,r,'edgecolor','none') %0:.1:1
caxis([0 1])
contour(R,I,r,0:.1:1.1) %0:.1:1
hold on

%% Exact roots
% mob=@(x)(x+1)./(x-1);
root=TransRoot(Dc,6,x0,Rad);

plot(root,'ro','markersize',6,'markerfacecolor','r')
plot([0 0],[-Rad Rad],'b','linewidth',1.5)
plot([-Rad Rad],[0 0],'k--','linewidth',1)
xlim([-Rad Rad])
ylim([-Rad Rad])

%%

lambda=1i*(-Rad:Rad/1000:Rad);
[~,Ind]=min(SpectralRadius(W,lambda,perturbation));
            if Ind==1
                Ind=2;
            elseif Ind==length(lambda)
                Ind=length(lambda)-1;
                
            end
            lref=1i*(imag(lambda(1,Ind-1)):imag(lambda(1,Ind+1)-lambda(1,Ind-1))/1000:imag(lambda(1,Ind+1)));
            
SR=SpectralRadius(W,lref,perturbation);
r0=min(SR)

contour(R,I,r,[r0 r0],'k') %0:.1:

% figure(1000);plot(abs(lambda),SR);