function Eq = PDmax_fval(ax,par)
%%
%
NumEq=1;
NumVar=size(ax,2);

Eq=zeros(NumEq,NumVar);

for kax=1:NumVar
    
    om=ax(1,kax);
    
    b=3*9.81/2/par.L;
    tau=par.tau;
    et=par.et;
    
    Eq(1,kax)=b.^(1/2).*(b+om.^2).*(b.^(1/2).*om.*(3+(-2).*cos(et.*om.* ...
  tau)+4.*cosh(b.^(1/2).*(1+et).*tau).*sin((1/2).*et.*om.*tau) ...
  .*sin((1/2).*(2+et).*om.*tau))+2.*(b+(-1).*om.^2).*cos((1/2) ...
  .*(2+et).*om.*tau).*sin((1/2).*et.*om.*tau).*sinh(b.^(1/2).* ...
  (1+et).*tau)).^(-1).*(om+2.*sin((1/2).*et.*om.*tau).*(om.* ...
  cosh(b.^(1/2).*(1+et).*tau).*sin((1/2).*(2+et).*om.*tau)+ ...
  b.^(1/2).*cos((1/2).*(2+et).*om.*tau).*sinh(b.^(1/2).*(1+et) ...
  .*tau)))-b;
end
end