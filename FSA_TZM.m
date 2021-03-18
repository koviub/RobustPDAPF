function Eq = FSA_TZM(ax,par)
%%
%
NumEq=1;
NumVar=size(ax,2);

Eq=zeros(NumEq,NumVar);

for kax=1:NumVar
    
    L=ax(2,kax);
    tau=ax(1,kax);
%     et=ax(3,kax);
    
    b=3/2*9.81/L;
    et=par.et;
    p=par.pp(0,b,tau,et);
    d=par.dd(0,b,tau,et);
    
    Eq(1,kax)=par.Dc2(0,b,tau,p,d,et);
end
end

