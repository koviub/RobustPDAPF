clc
%% invetistigated domain
x0=0+1i*0;                          % origin on the complex plane
if strcomp(uncertainty,'delay')
    Rad=6;
else
    Rad=30;
end                            % radius around x0
lambda=1i*(0:Rad/1000:Rad);

%% Properties
uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY'
structure    = 'unSTRUCTURED'; % 'STRUCTURED'/'UNSTRUCTURED'
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();

step0=30;step1=20;
Radii=zeros(step0);Area=zeros(step0); % stability radius, size of stable region [unit^2/s^3] 

% Parameters for domain (tau,L)
Ls=0.01;Le=10;ts=.05;te=1;
tt=ts:(te-ts)/(step0-1):te;
ll=Ls:(Le-Ls)/(step0-1):Le;
[tk,lk]=meshgrid(tt,ll);
for s=1:length(lk)
lk(s,:)=lk(s,:)+3/4*9.81*tk(s,:).^2;
end

fprintf('progres: %.1f percent \n',0);
for m=1:step0
    for k=1:step0
        
        par.L=lk(k,m);
        par.tau=tk(k,m);
        
        [ps,pe,ds,de]=PDRegion(par);
        
        stabpoint=0;
        r=zeros(step1);
        pp=ps:(pe-ps)/(step1-1):pe;
        dd=ds:(de-ds)/(step1-1):de;
        [pk,dk]=meshgrid(pp,dd);
        for j=1:step1
            for i=1:step1
                fprintf('*');
                par.p=pk(i,j);
                par.d=dk(i,j);
                
                [W,Mx,Dc,w,G]=ModelDefinition(par,uncertainty,structure);
                
                val=max(abs(eig(G)));
                
                if val<=1
                    r(i,j)=min(SpectralRadius(W,lambda,perturbation));
                    stabpoint=stabpoint+1;
                else
                    r(i,j)=0;
                end
                
            end
            fprintf('\n')
        end
        clc
        Area(k,m)=stabpoint/step1^2*(pe-ps)*(de-ds); 
        Radii(k,m)=max(max(r));
        fprintf('progres %.1f percent \n',((m-1)*step0+k)/step0^2*100)
        
    end
end

figure(10)
hold on;grid on;
plot(tt,3/(4*(1+par.a))*9.81*tt.^2,'b')
surface(tk,lk,Radii);
ylim([0 10])
caxis([0 1])
xlabel('$\tau$','interpreter','latex')
ylabel('$r^{\delta}$','interpreter','latex')

figure(11)
hold on;grid on;
surface(tk,lk,log(Area))
plot(tt,3/4*9.81.*tt.^2,'r')
ylim([0 10])
caxis([-13 8])
xlabel('$\tau$','interpreter','latex')
ylabel('Area','interpreter','latex')

res.Radii=Radii;
res.Area=Area;
res.tau=tk;
res.Length=lk;

save filname.mat res
