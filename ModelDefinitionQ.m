function [W,Mx,Dc,w,w1,w2,G]=ModelDefinitionQ(par,str,str2)
%% system parameters
L=par.L;
g=9.81;
b=3/2*g/L;

tau=par.tau;
p=par.p;
d=par.d;
a=par.a;

q=par.q;

% h=par.h;
% rr=floor(tau/h);
rr=par.r;
h=tau/rr;

% state-space representation
if q~=0
    A=[0,1;...
        b/q,(q-1)/q];
    B=[0;...
        1];
    C1=[-p/q,-d/q];
    C2=[0,0];
else
    A=b;
    B=1;
    C1=-p;
    C2=-d;
    
end

% Semi-discretization
P1=expm(A*h);
R1=(P1-eye(length(P1)))/A*B;
P2=A*P1;
R2=(A*R1+B);

N=length(P1);
G=zeros(2*N+rr);
for i=2:2*N+rr
    G(i,i-1)=1;
end
G(1:N,1:N)=P1;
G(1:N,end)=R1;
G(N+1:2*N,1:N)=P2;
G(N+1:2*N,end)=R2;
G(2*N+1,1:N)=C1;
G(2*N+1,N+1:2*N)=C2;

% Characteristic matrix and equation
Mx=@(x)x*eye(length(A))-A-B*C1*exp(-x*tau)-x*B*C2*exp(-x*tau);
Dc=@(x)x.^2*q+x.*(1-q+d*exp(-x*tau))+(p*exp(-x*tau)-b);
%det(Mx(lambda));

%% Problem definition
switch lower(str)
    case 'system'
        w=@(x)eye(length(A));
        w1=@(x)-b;
        w2=@(x)0;
    case 'control'
        Mx=@(x)transpose(Mx(x));
        w=@(x)transpose((C1*tau+x*C2)*exp(-x));
        w1=@(x)exp(-x*tau)*[p d*x];
        w2=@(x)0;
        
    case 'delay'
        % -----------
        % in progress
        % ???????????
        w=@(x)A+B*C1*exp(-x);
        w1=@(x)x.*d*exp(-(x))+2*(p*exp(-(x))-b);
        w2=@(x)p*exp(-(x))-b;
        
end

%         Dc=@(x)Dc(tau*x);
%         w=@(x)w(tau*x);
%         w1=@(x)w1(tau*x);
%         w2=@(x)w2(tau*x);

switch lower(str2)
    case 'structured'
        W=@(x)w1(x)./Dc(x);
        
    case 'unstructured'
        W=@(x)Mx(x)\w(x);
end
end