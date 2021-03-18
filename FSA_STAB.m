function Eq = FSA_STAB(ax,par)
%%
%
NumEq=2;
NumVar=size(ax,2);

Eq=zeros(NumEq,NumVar);

for kax=1:NumVar
    
    p=ax(1,kax);
    d=ax(2,kax);
    om=ax(3,kax);

    et=par.et;
    L=par.L;
    tau=par.tau;
    b=3/2*9.81/L;
    
    lam=1i*om;

    
%     a1=-1/2*sqrt(b)*exp(-sqrt(b)*(tau+et*tau))+1/2*sqrt(b)*exp(sqrt(b)*(tau+et*tau));
%     a2=1/2*exp(-sqrt(b)*(tau+et*tau))+1/2*exp(sqrt(b)*(tau+et*tau));
%     a3=(exp(sqrt(b)*(tau+et*tau))-exp(-sqrt(b)*(tau+et*tau)))/(2*sqrt(b));
%     a4=1/2*((-1+exp(tau+et*tau)*(sqrt(b)-lam))/(sqrt(b)-lam)-(-1+exp(tau+et*tau)*(-sqrt(b)-lam))/(sqrt(b)+lam));
%     a5=1/2*(-(1-exp((tau+et*tau)*(-sqrt(b)-lam)))/(sqrt(b)*(sqrt(b)+lam))+(-1+exp((tau+et*tau)*(sqrt(b)-lam)))/(b-sqrt(b)*lam));
%     T=[[lam,-1,0];...
%         [-b,lam,-exp(-lam*tau)];...
%         [a1*d+a2*p,a2*d+a3*p,1-a4*d-a5*p]];
    Dcar=CharEq(lam,b,tau,p,d,et);
    
    Eq(1,kax)=real(Dcar);
    Eq(2,kax)=imag(Dcar);
end
end

function DC=CharEq(lam,b,tau,p,d,et)

DC=-b + lam*(d + lam) + p + (2*(-(sqrt(b)*(d*lam + p)*cosh(sqrt(b)*(1 + et)*tau)) - (b*d + lam*p)*sinh(sqrt(b)*(1 + et)*tau))*sinh((lam*(tau - (1 + et)*tau))/2))/(sqrt(b)*exp((lam*(tau + (1 + et)*tau))/2));
    
end