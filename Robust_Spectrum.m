clear all
% clc

uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY/'ALL''
structure    = 'STRUCTURED'; % 'STRUCTURED'/'UNSTRUCTURED'
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();

%% modell

[W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,structure);
figure();
for ii=1:4
%% invetistigated domain
x0=0+1i*(ii-1)*60;                          % origin on the complex plane
Rad=30;                         % radius around x0

re=-Rad:Rad/200:Rad;
im=-Rad:Rad/200:Rad;
[R,I]=meshgrid(re,im);
lambda=R+1i*I;

%% Exact roots
% mob=@(x)(x+1)./(x-1);
% figure();

root=TransRoot(Dc,8,x0,Rad);


hold on
plot(root,'ro','markersize',6,'markerfacecolor','r')
plot([0 0],[-Rad Rad],'b','linewidth',1.5)
plot([-Rad Rad],[0 0],'k--','linewidth',1)
% xlim([-Rad Rad])
% ylim([-Rad Rad])

%% Robustness

r=SpectralRadius(W,lambda,perturbation);


% surface(R,I,r,'edgecolor','none') %0:.1:1
caxis([0 1])
contour(R,I,r) %0:.1:1
drawnow();
end

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