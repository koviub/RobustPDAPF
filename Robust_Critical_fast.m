clc
%% invetistigated domain
x0=0+1i*0;                          % origin on the complex plane
% if strcomp(uncertainty,'delay')
%     Rad=6;
% else
%     Rad=30;
% end                             % radius around x0
Rad=30;
lambda=1i*(0:Rad/200:Rad);

%% Properties
uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY'
structure    = 'STRUCTURED'; % 'STRUCTURED'/'UNSTRUCTURED'
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();

step0=50;step1=50;
Radii=zeros(step0);Area=zeros(step0); % stability radius, size of stable region [unit^2/s^3]

% Parameters for domain (tau,L)
Ls=0.01;Le=10;ts=.05;te=1;
tt=ts:(te-ts)/(step0-1):te;
ll=Ls:(Le-Ls)/(step0-1):Le;
[tk,lk]=meshgrid(tt,ll);
% for s=1:length(lk)
%     lk(s,:)=lk(s,:)+3/4*9.81*tk(s,:).^2;
% end

fprintf('progres: %.1f percent \n',0);
for m=1:step0
    for k=1:step0
        par.L=lk(k,m);
        par.tau=tk(k,m);
        
        [ps,pe,ds,de]=PDRegion(par);
        
        val=zeros(step1);r=zeros(step1);
        pp=ps:(pe-ps)/(step1-1):pe;
        dd=ds:(de-ds)/(step1-1):de;
        [pk,dk]=meshgrid(pp,dd);
        
        i0=ceil(step1/2);inext=0;
        j0=ceil(step1/2);jnext=0; % system 3, control 4
        i1=[0 1 1 1 0 -1 -1 -1 0 1 2 2 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1]; % 0 1 2 2 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1
        j1=[1 1 0 -1 -1 -1 0 1 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1 0 1 2 2]; % 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1 0 1 2 2
        
        % first point
        par.p=pk(i0,j0);
        par.d=dk(i0,j0);
        
        [W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,structure);
        
        val(i0,j0)=max(abs(eig(G)));
        
        if val(i0,j0)<=1
            %             r(i0,j0)=min(SpectralRadius(W,lambda,perturbation));
            [~,Ind]=min(SpectralRadius(W,lambda,perturbation));
            if Ind==1
                Ind=2;
            elseif Ind==length(lambda)
                Ind=length(lambda)-1;
                
            end
            lref=1i*(imag(lambda(1,Ind-1)):imag(lambda(1,Ind+1)-lambda(1,Ind-1))/1000:imag(lambda(1,Ind+1)));
            r(i0,j0)=min(SpectralRadius(W,lref,perturbation));
            
        else
            r(i0,j0)=0;
        end
        
        % find maximum
        while true
            for f=1:length(i1)
                indi=i0+i1(f);
                indj=j0+j1(f);
                
                if indi<1
                    indi=1;
                elseif indi>step1
                    indi=step1;
                end
                
                if indj<1
                    indj=1;
                elseif indj>step1
                    indj=step1;
                end
                
                if r(indi,indj)==0
                    fprintf('*');
                    par.p=pk(indi,indj);
                    par.d=dk(indi,indj);
                    
                    [W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,structure);
                    
                    val(indi,indj)=max(abs(eig(G)));
                    if val(indi,indj)<=1
                        %                         r(indi,indj)=min(SpectralRadius(W,lambda,perturbation));
                        [~,Ind]=min(SpectralRadius(W,lambda,perturbation));
                        if Ind==1
                            Ind=2;
                        elseif Ind==length(lambda)
                            Ind=length(lambda)-1;
                            
                        end
                        lref=1i*(imag(lambda(1,Ind-1)):imag(lambda(1,Ind+1)-lambda(1,Ind-1))/1000:imag(lambda(1,Ind+1)));
                        r(indi,indj)=min(SpectralRadius(W,lref,perturbation));
                        
                    else
                        r(indi,indj)=0;
                    end
                end
                
                if r(i0,j0)<r(indi,indj)
                    inext=indi;jnext=indj;
                end
            end
            fprintf('\n');
            if inext==i0&&jnext==j0
                rmax=r(i0,j0);
                fprintf('found max');
                break;
                
            end
            
            if inext<1
                inext=1;
            elseif inext>step1
                inext=step1;
            end
            
            if jnext<1
                jnext=1;
            elseif jnext>step1
                jnext=step1;
            end
            i0=inext;j0=jnext;
        end
        clc
        Radii(k,m)=rmax;
        fprintf('progres %.1f percent \n',((m-1)*step0+k)/step0^2*100);
    end
end

figure(11)
hold on;grid on;
% plot(tt,3/(4*(1+par.a))*9.81*tt.^2,'b')
% plot(tt,3./(4*(par.q+(1-par.q).*tt))*9.81*tt.^2,'b')
surface(tk,lk,Radii);
ylim([0 10])
caxis([0 1])
xlabel('$\tau$','interpreter','latex')
ylabel('$r^{\delta}$','interpreter','latex')

res.Radii=Radii;
res.Area=Area;
res.tau=tk;
res.Length=lk;

save(lower(sprintf('%s_%s_%s_PDA_07.mat',uncertainty,perturbation,structure)),'res')
% save filname.mat res
