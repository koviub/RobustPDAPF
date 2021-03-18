function [W,Mx,Dc,w,w1,w2,G]=ModelDefinition1(par,str,str2)
%% system parameters
L=par.L;
b=3/2*9.81/L;

tau=par.tau;
p=par.p;
d=par.d;

et=par.et;
tauh=tau*(1+et);

% h=par.h;
% rr=floor(tau/h);
rr=par.r;
h=tau/rr;

rrh=floor(tauh/h);
% state-space representation
A=[0,1;...
    b,0];
B=[0;...
    1];
C1=[-p,-d];

% Semi-discretization
P1=expm(A*h);
R1=(P1-eye(length(P1)))/A*B;
F=C1*expm(A*tauh);
Q=@(j)C1*expm(A*j*h)*B*h;

N=length(P1);
rs=max(rr,rrh);
G=zeros(N+rs);
for i=2:N+rs
    G(i,i-1)=1;
end
G(1:N,1:N)=P1;
G(1:N,N+rr)=R1;
G(N+1,1:N)=F;

for ii=1:rrh
    G(N+1,N+ii)=Q(ii);
end

M1=@(x)x*eye(N)-A;
M2=@(x)-B*exp(-x*tau);
M3=F;
M4=@(x)1+(1/2).*d.*(((-1)+exp(1).^((1+et).*tau.*(b.^(1/2)+(-1).*x)) ...
    ).*(b.^(1/2)+(-1).*x).^(-1)+(-1).*((-1)+exp(1).^((-1).*(1+ ...
    et).*tau.*(b.^(1/2)+x))).*(b.^(1/2)+x).^(-1))+(1/2).*p.*(( ...
    -1).*b.^(-1/2).*(1+(-1).*exp(1).^((-1).*(1+et).*tau.*(b.^( ...
    1/2)+x))).*(b.^(1/2)+x).^(-1)+((-1)+exp(1).^((1+et).*tau.*( ...
    b.^(1/2)+(-1).*x))).*(b+(-1).*b.^(1/2).*x).^(-1));
% Characteristic matrix and equation
Mx=@(x)[M1(x) M2(x);M3 M4(x)];
Dc=@(x)(-1).*b+x.*(d+x)+p+2.*b.^(-1/2).*exp(1).^((-1/2).*(2+et) ...
    .*x.*tau).*(b.^(1/2).*(d.*x+p).*cosh(b.^(1/2).*(1+et).* ...
    tau)+(b.*d+x.*p).*sinh(b.^(1/2).*(1+et).*tau)).*sinh((1/2) ...
    .*et.*x.*tau);
%det(Mx(lambda));

%% Problem definition
switch lower(str)
    case 'system'
        w=@(x)eye(length(A));
        w1=@(x)-b*tau^2;
        %         w1=@(x)[x^2, -b*tau^2];
        w2=@(x)0;
    case 'control'
        Mx=@(x)transpose(Mx(x));
        w=@(x)transpose((C1*tau+x*C2)*exp(-x));
        w1=@(x)[1+2.*exp(1).^((-1/2).*(2+et).*tau.*x).*cosh(b.^(1/2).*(1+et) ...
            .*tau).*sinh((1/2).*et.*tau.*x)+2.*b.^(-1/2).*exp(1).^(( ...
            -1/2).*(2+et).*tau.*x).*x.*sinh(b.^(1/2).*(1+et).*tau).* ...
            sinh((1/2).*et.*tau.*x);...
            x+2.*exp(1).^((-1/2).*(2+et).*tau.*x).*x.*cosh(b.^(1/2).*(1+ ...
            et).*tau).*sinh((1/2).*et.*tau.*x)+2.*b.^(1/2).*exp(1).^(( ...
            -1/2).*(2+et).*tau.*x).*sinh(b.^(1/2).*(1+et).*tau).*sinh(( ...
            1/2).*et.*tau.*x)];

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