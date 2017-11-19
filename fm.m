function [sp,sn,XPnoise,XNnoise,ET]= fm(Xp,Xn,Parameter)
%  Author: Bin-Bin Gao 
%  Email: csgaobb@gmail.com
%  July 5, 2017

% check correct number of arguments

if ( nargin>3||nargin<3) 
    help fm
else
    if (nargin<3) Parameter.ker='linear'; end
    Parameter.u = 0.1;
    u=Parameter.u;
    eplison=1e-10;
    [rxp,~]=size(Xp);
    [rxn,~]=size(Xn);
    
    switch Parameter.ker
        case 'linear'
            kfun = @linear_kernel;kfunargs ={};
        case 'quadratic'
            kfun = @quadratic_kernel;kfunargs={};
        case 'radial'
            p1=Parameter.p1;
            kfun = @rbf_kernel;kfunargs = {p1};
        case 'rbf'
            p1=Parameter.p1;
            kfun = @rbf_kernel;kfunargs = {p1};
        case 'polynomial'
            p1=Parameter.p1;
            kfun = @poly_kernel;kfunargs = {p1};
        case 'mlp'
            p1=Parameter.p1;
            p2=Parameter.p2;
            kfun = @mlp_kernel;kfunargs = {p1, p2};
    end
    switch lower(Parameter.ker)
        case 'linear'
            stt = cputime;
            Xp_cen=mean(Xp);%
            Xn_cen=mean(Xn);%
            % radiusxp
            radiusxp=sum((repmat(Xp_cen,rxp,1)-Xp).^2,2);
            radiusxpxn=sum((repmat(Xn_cen,rxp,1)-Xp).^2,2);
          
            radiusmaxxp=max(radiusxp);
            % radiusxn
            radiusxn=sum((repmat(Xn_cen,rxn,1)-Xn).^2,2);%||Xi--Xcen-||^2
            radiusxnxp=sum((repmat(Xp_cen,rxn,1)-Xn).^2,2);%||xi--Xcen+||^2
           
            radiusmaxxn=max(radiusxn);
           
            sp=zeros(rxp,1);
            XPnoise=find(radiusxp>=radiusxpxn);
            XPnormal=find(radiusxp<radiusxpxn);
            sp(XPnormal,1)=(1-u).*(1-sqrt(abs(radiusxp(XPnormal,1))./(radiusmaxxp+eplison)));
            sp(XPnoise,1)=u.*(1-sqrt(abs(radiusxp(XPnoise,1))./(radiusmaxxp+eplison)));
            
           
            sn=zeros(rxn,1);
            XNnoise=find(radiusxn>=radiusxnxp);
            XNnormal=find(radiusxn<radiusxnxp);
            sn(XNnormal,1)=(1-u)*(1-sqrt(abs(radiusxn(XNnormal,1))./(radiusmaxxn+eplison)));
            sn(XNnoise,1)=u.*(1-sqrt(abs(radiusxn(XNnoise,1))./(radiusmaxxn+eplison)));
            
            ET=cputime - stt;            
        case 'rbf'
            sttt = cputime;
            %K(xi+,xj+)
            kerxp1 = feval(kfun,Xp,Xp,kfunargs{:});%K(X+,X+)
            % K(xi-,xj-)
            kerxp2 = feval(kfun,Xn,Xn,kfunargs{:});%K(X-,X-)
            % K(xi+,xj-)
            kerxp3 = feval(kfun,Xp,Xn,kfunargs{:});%K(X+,X-)
            % K(xi-,xj+)==K(xi+,xj-)^T
            radiusxp=1-2*mean(kerxp1,2)+mean(mean(kerxp1));%||(xi+)-cen+||^2
            radiusxpxn=1-2*mean(kerxp3,2)+mean(mean(kerxp2));%  ||(xi+)-cen-||^2
            radiusmaxxp=max(radiusxp);
            
            radiusxn=1-2*mean(kerxp2,2)+mean(mean(kerxp2));   %  ||(xi-)-cen-||^2
            radiusxnxp=1-2*transpose(mean(kerxp3,1))+mean(mean(kerxp1)); %  ||(xi-)-cen+||^2
            radiusmaxxn=max(radiusxn);
            sp=zeros(rxp,1);
            XPnoise=find(radiusxp>=radiusxpxn);
            XPnormal=find(radiusxp<radiusxpxn);
            sp(XPnormal,1)=(1-u).*(1-sqrt(abs(radiusxp(XPnormal,1))./(radiusmaxxp+eplison)));
            sp(XPnoise,1)=u.*(1-sqrt(abs(radiusxp(XPnoise,1))./(radiusmaxxp+eplison)));           
            sn=zeros(rxn,1);
            XNnoise=find(radiusxn>=radiusxnxp);
            XNnormal=find(radiusxn<radiusxnxp);
            sn(XNnormal,1)=(1-u)*(1-sqrt(abs(radiusxn(XNnormal,1))./(radiusmaxxn+eplison)));
            sn(XNnoise,1)=u.*(1-sqrt(abs(radiusxn(XNnoise,1))./(radiusmaxxn+eplison)));
            ET=cputime - sttt;
    end
      sp=mapminmax(sp',eps,1)';
      sn=mapminmax(sn',eps,1)';
end

