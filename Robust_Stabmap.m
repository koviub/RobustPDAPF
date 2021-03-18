clc
%% invetistigated domain
x0=0+1i*0;                          % origin on the complex plane
% if strcomp(uncertainty,'delay')
%     Rad=6;
% else
% end                            % radius around x0
Rad=30;
lambda=1i*(0.00:Rad/5000:Rad);

%% Properties
uncertainty  = 'control';    % 'SYSTEM'/'CONTROL'/'DELAY'
structure    = 'STRUCTURED'; % 'STRUCTURED'/'UNSTRUCTURED'
perturbation = 'real';       % 'REAL'/'COMPLEX'

par=Parameters();
% figure();
[ps,pe,ds,de]=PDRegion(par,12);

step=50;
val=zeros(step);r=zeros(step);
pp=ps:(pe-ps)/(step-1):pe;
dd=ds:(de-ds)/(step-1):de;
[pk,dk]=meshgrid(pp,dd);
sf=surface(pk,dk,r,'edgecolor','none');
for j=1:step
    fprintf('%d\n',j)
    for i=1:step
        
        par.p=pk(i,j);
        par.d=dk(i,j);
        
        [W,Mx,Dc,w,w1,w2,G]=ModelDefinition(par,uncertainty,structure);
        
        val(i,j)=max(abs(eig(G)));
        if val(i,j)<=1
            
            [~,Ind]=min(SpectralRadius(W,lambda,perturbation));
            if Ind==1
                Ind=2;
            elseif Ind==length(lambda)
                Ind=length(lambda)-1;
                
            end
            lref=1i*(imag(lambda(1,Ind-1)):imag(lambda(1,Ind+1)-lambda(1,Ind-1))/1000:imag(lambda(1,Ind+1)));
            [r(i,j),~]=min(SpectralRadius(W,lref,perturbation));
            %
            %             SR=SpectralRadius(W,lambda,perturbation);
            % %             figure(1000);plot(abs(lambda),SR);
            %             [r(i,j),~]=min(SR);
            % %             fprintf('%d',ii)
        else
            r(i,j)=0;
        end
        
        set(sf,'ZData',r,'CData',r);
        drawnow();
    end
end
set(sf,'Visible','off')
contour(pk,dk,r,0:.1:1);
