function Eq = PDmax_fval_q(ax,par)
%%
%
NumEq=1;
NumVar=size(ax,2);

Eq=zeros(NumEq,NumVar);

for kax=1:NumVar
    
    om=ax(1,kax);
    
    II=par.In;
    GG=par.Gr;
    b=GG/II;
    tau=par.tau;
    q=par.q;
    
    Eq(1,kax)=(b+q.*om.^2).*cos(om.*tau)+(1-q).*om.*sin(om.*tau);
end
end