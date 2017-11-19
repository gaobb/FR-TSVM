function  [frtsvm_struct] = frtsvmtrain(Traindata,Trainlabel,Parameter)
%  Author: Bin-Bin Gao 
%  Email:csgaobb@gmail.com
%  July 5, 2017

% check correct number of arguments
if ( nargin>3||nargin<3) 
    help  frtsvmtrain
end

ker=Parameter.ker;
CC=Parameter.CC;
CR=Parameter.CR;
% Parameter.autoScale=0;
showplots=Parameter.showplots;
autoScale=0;

% 1 
tic;
[groupIndex, groupString] = grp2idx(Trainlabel);
groupIndex = 1 - (2* (groupIndex-1));
scaleData = [];
if autoScale
    scaleData.shift = - mean(Traindata);
    stdVals = std(Traindata);
    scaleData.scaleFactor = 1./stdVals;
    % leave zero-variance data unscaled:
    scaleData.scaleFactor(~isfinite(scaleData.scaleFactor)) = 1;
    % shift and scale columns of data matrix:
    for k = 1:size(Traindata, 2)
        scTraindata(:,k) = scaleData.scaleFactor(k) * ...
            (Traindata(:,k) +  scaleData.shift(k));
    end
else
    scTraindata= Traindata;
end

Xp=scTraindata(groupIndex==1,:);
Lp=Trainlabel(groupIndex==1);
Xn=scTraindata(groupIndex==-1,:);
Ln=Trainlabel(groupIndex==-1);
X=[Xp;Xn];
L=[Lp;Ln];
[sp,sn,NXpv,NXnv]=fm(Xp,Xn,Parameter);

lp=sum(groupIndex==1);
ln=sum(groupIndex==-1);
switch ker
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
switch ker
    case 'linear'
        Kpx=Xp;Knx=Xn;
    case 'rbf'
        Kpx = feval(kfun,Xp,X,kfunargs{:});%K(X+,X)
        Knx = feval(kfun,Xn,X,kfunargs{:});%K(X-,X)
end
S=[Kpx ones(lp,1)];R=[Knx ones(ln,1)];

CC1=CC*sn;
CC2=CC*sp;

fprintf('Optimising ...\n');
switch  Parameter.algorithm
    case  'cd'
        [alpha, vp] = cd(S,R,CR,CC1);
        [beta , vn] = cd(R,S,CR,CC2);
        vn=-vn;
    case  'qp'
        QR=(S'*S+CR*eye(size(S'*S)))\R';
        RQR=R*QR;
        RQR=(RQR+RQR')/2;
        
        QS=(R'*R+CR*eye(size(R'*R)))\S';
        SQS=S*QS;
        SQS=(SQS+SQS')/2;

        [alpha,~,~]=qp(RQR,-ones(ln,1),[],[],zeros(ln,1),CC1,ones(ln,1));
        [beta,~,~] =qp(SQS,-ones(lp,1),[],[],zeros(lp,1),CC2,ones(lp,1));
        
        vp=-QR*alpha;
        vn=QS*beta;
end
Time= toc;
frtsvm_struct.scaleData=scaleData;

frtsvm_struct.X = X;
frtsvm_struct.L = L;
frtsvm_struct.sp = sp;
frtsvm_struct.sn = sn;

frtsvm_struct.alpha = alpha;
frtsvm_struct.beta  = beta;
frtsvm_struct.vp = vp;
frtsvm_struct.vn = vn;

frtsvm_struct.KernelFunction = kfun;
frtsvm_struct.KernelFunctionArgs = kfunargs;
frtsvm_struct.Parameter = Parameter;
frtsvm_struct.groupString=groupString;
frtsvm_struct.time=Time;

frtsvm_struct.NXpv=NXpv;
frtsvm_struct.NXnv=NXnv;
frtsvm_struct.nv=length(NXpv)+length(NXnv);
if  showplots
    frtsvmplot(frtsvm_struct,Traindata,Trainlabel);
end   
end






