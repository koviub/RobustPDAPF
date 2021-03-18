clear all
clc
%% invetistigated domain
x0=0+1i*0;                          % origin on the complex plane
% if strcomp(uncertainty,'delay')
%     Rad=6;
% else
%     Rad=30;
% end      % radius around x0
Rad=10000;
lambda=1i*(0:Rad/50000:Rad);

%% Properties
uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY'
structure    = 'STRUCTURED'; % 'STRUCTURED'/'UNSTRUCTURED'
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();
figure
[ps,pe,ds,de]=PDRegion(par,105);

step=50;
val=zeros(step);r=zeros(step);
pp=ps:(pe-ps)/(step-1):pe;
dd=ds:(de-ds)/(step-1):de;
[pk,dk]=meshgrid(pp,dd);
sf=surface(pk,dk,r);

i0=ceil(step/2);inext=i0-1;
j0=ceil(step/2);jnext=j0-1;% control 4, system 3
i1=[0 1 1 1 0 -1 -1 -1 0 1 2 2 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1]; % 0 1 2 2 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1
j1=[1 1 0 -1 -1 -1 0 1 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1 0 1 2 2]; % 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1 0 1 2 2

% first point
par.p=pk(i0,j0);
par.d=dk(i0,j0);

[W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,structure);

val(i0,j0)=max(abs(eig(G)));
if val(i0,j0)<=1
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
set(sf,'ZData',r,'CData',r,'edgecolor','none');
drawnow();

k=0;
while true% found maximum
    
    for f=1:length(i1)
        indi=i0+i1(f);
        indj=j0+j1(f);
        
        if indi<1
            indi=1;
        elseif indi>step
            indi=step;
        end
        
        if indj<1
            indj=1;
        elseif indj>step
            indj=step;
        end
        
        if r(indi,indj)==0
            
            k=k+1;
            fprintf('box counted: %d \n',k);
            par.p=pk(indi,indj);
            par.d=dk(indi,indj);
            
            [W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,structure);
            
            val(indi,indj)=max(abs(eig(G)));
            if val(indi,indj)<=1
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
            
            set(sf,'ZData',r,'CData',r);
            drawnow();
        end
        
        if r(i0,j0)<r(indi,indj)
            inext=indi;jnext=indj;
            %             fprintf('%f,%f',r(i0,j0),r(indi,indj))
        end
        
    end
    
    if inext==i0&&jnext==j0
        fprintf('Maximum found! r=%f \n',r(i0,j0));
        break;
    end
    i0=inext;j0=jnext;
    %     k=k+1;
    %     fprintf('iteration %d: %d,%d \n',k,i0,j0);
end
