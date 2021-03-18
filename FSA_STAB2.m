function Eq = FSA_STAB2(ax,par)
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

    Dcar=(-((sqrt(b) - lam)*(sqrt(b)*d - p)) + 2*sqrt(b)*exp((1 + et)*(sqrt(b) + lam)*tau)*(b - lam*(d + lam) - p) + exp(2*sqrt(b)*(1 + et)*tau)*(sqrt(b) + lam)*(sqrt(b)*d + p))/(2*sqrt(b)*exp((1 + et)*(sqrt(b) + lam)*tau)*(b - lam^2));
    Eq(1,kax)=real(Dcar);
    Eq(2,kax)=imag(Dcar);
end
end