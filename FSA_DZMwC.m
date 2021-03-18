function Eq = FSA_DZMwC(ax,par)
%%
%
NumEq=2;
NumVar=size(ax,2);

Eq=zeros(NumEq,NumVar);

for kax=1:NumVar
    
    tau=ax(1,kax);
    L=ax(2,kax);
    om=ax(3,kax);
    
    b=3/2*9.81/L;
    et=par.et;
    p=par.pp(0,b,tau,et);
    d=par.dd(0,b,tau,et);
    
    Eq(1,kax)=real(par.Dc0(1i*om,b,tau,p,d,et));
    Eq(2,kax)=imag(par.Dc0(1i*om,b,tau,p,d,et));
end
end

