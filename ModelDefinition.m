function [W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,str,str2)
%% system parameters
L=par.L;
b=3/2*9.81/L;

tau=par.tau;
p=par.p;
d=par.d;
a=par.a;

% h=par.h;
% rr=floor(tau/h);
rr=par.r;
h=tau/rr;

% state-space representation
A=[0,1;...
    b,0];
B=[0;...
    1];
C1=[-p,-d];
C2=[0,-a];

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
Dc=@(x)x.^2.*(1+a*exp(-x.*tau))+x.*d*exp(-x.*tau)+p*exp(-x.*tau)-b;
%det(Mx(lambda));

%% Problem definition
switch lower(str)
    case 'system'
        Wu=@(x)inv(Mx(x));
        w=@(x)b;
        Ws=@(x)w(x)./Dc(x);
        w1=@(x)0;
        w2=@(x)0;
    case 'control'
        Wu=@(x)(C1+x*C2)/Mx(x)*exp(-x*tau);
        w=@(x)exp(-x*tau)*...
            [p;...
            d*x;...
            a*x^2];
        Ws=@(x)w(x)./Dc(x);
        w1=@(x)0;
        w2=@(x)0;
    case 'all'
        Wu=@(x)(C1+x*C2)/Mx(x)*exp(-x*tau);
        w=@(x)exp(-x*tau)*...
            [0;p;...
            d*x;...
            a*x^2]+[-b;0;0;0];
        Ws=@(x)w(x)./Dc(x);
        w1=@(x)0;
        w2=@(x)0;
    case 'delay'
        % -----------
        % in progress
        % ???????????
        
        Mx=@(x)x*eye(length(A))-tau*A-B*(C1*tau-x*C2)*exp(-x);
        w=@(x)A+B*C1*exp(-x);
        Dc=@(x)x.^2*(1+a*exp(-x))+x.*d*exp(-x).*tau+(p*exp(-x)-b).*tau^2;
        w1=@(x)x.*d*exp(-(x))+2*(p*exp(-(x))-b);
        w2=@(x)p*exp(-(x))-b;
        
        %         Dc=@(x)Dc(tau*x);
        %         w1=@(x)w1(tau*x);
        %         w2=@(x)w2(tau*x);
        Wu=@(x)Mx(tau*x)\w(tau*x);
        Ws=@(x)(w1(x)+w2(x)*0)./Dc(x);
end

switch lower(str2)
    case 'structured'
        W=@(x)Ws(x);
    case 'unstructured'
        W=@(x)Wu(x);
end
end