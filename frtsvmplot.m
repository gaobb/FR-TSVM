function frtsvmplot(frtsvm_struct,Traindata,Trainlabel)
% Function:  visualizing 
% Input:      
% frtsvm_struct,Traindata,Trainlabel

% Output:    
%  Author: Bin-BinGaa (csgaobb@gmail.com)
% Created on 2016.10.10
% Last modified on 2017.07.16
if (nargin <1|| nargin > 3) % check correct number of arguments
    help fftsvmplot
else

% if ~isempty(frtsvm_struct.scaleData)
%     scaleData=frtsvm_struct.scaleData;
%     for c = 1:size(Traindata, 2)
%         Traindata(:,c) = scaleData.scaleFactor(c) * ...
%             (Traindata(:,c) +  scaleData.shift(c));
%     end
% end    
    
    
%%
minX = min(Traindata(:, 1));
maxX = max(Traindata(:, 1));
minY = min(Traindata(:, 2));
maxY = max(Traindata(:, 2));

gridX = (maxX - minX) ./ 200;
gridY = (maxY - minY) ./ 200;

minX = minX - 10 * gridX;
maxX = maxX + 10 * gridX;
minY = minY - 10 * gridY;
maxY = maxY + 10 * gridY;

[bigX, bigY] = meshgrid(minX:gridX:maxX, minY:gridY:maxY);

%%

ntest=size(bigX, 1) * size(bigX, 2);
test_Traindata=[reshape(bigX, ntest, 1), reshape(bigY, ntest, 1)];
test_label = zeros(size(test_Traindata,1), 1);

[~,outclass,Z,f1,f2]= frtsvmclass(frtsvm_struct,test_Traindata);

bigC = reshape(outclass, size(bigX, 1), size(bigX, 2));
bigZ = reshape(Z, size(bigX, 1), size(bigX, 2));
bigf1 = reshape(f1, size(bigX, 1), size(bigX, 2));
bigf2 = reshape(f2, size(bigX, 1), size(bigX, 2));

alpha=frtsvm_struct.alpha;
beta =frtsvm_struct.beta;

[groupIndex, groupString] = grp2idx(Trainlabel);
groupIndex = 1 - (2* (groupIndex-1));
xp=Traindata(groupIndex==1,:);%�����ѵ����?
xn=Traindata(groupIndex==-1,:);%�����ѵ����frtsvm_struct.NXpv


ln=size(alpha,1);
lp=size(beta,1);

nsvIndex =  alpha > (frtsvm_struct.Parameter.CC*1e-4);%& alpha< (frtsvm_struct.Parameter.CC*frtsvm_struct.sn-sqrt(eps));
psvIndex =  beta > (frtsvm_struct.Parameter.CC*1e-4);%&  beta< (frtsvm_struct.Parameter.CC*frtsvm_struct.sp-sqrt(eps));

psv = xp(psvIndex,:);
nsv = xn(nsvIndex,:);
if ~size(psv,1)
    psvIndex =  beta > eps;
    psv = xp(psvIndex,:);
end
if ~size(nsv,1)
    nsvIndex =  alpha > eps;
    nsv = xn(nsvIndex,:);
end


figure
sp=mapminmax(frtsvm_struct.sp',0.1,1)';
[xq,yq] = meshgrid(linspace(min(xp(:, 1)), max(xp(:, 1)), 100), linspace(min(xp(:, 2)), max(xp(:, 2)), 100));
vq = griddata(xp(:, 1), xp(:, 2),sp,xq,yq,'v4');
% vq(find(vq<0))=0;
contourf(xq,yq,vq,100,'LineStyle','none')
colormap(jet);
hold on
plot(xp(:, 1), xp(:, 2),'r+','LineWidth',2);
colorbar
hold off
% line=['Fuzzy membership distribution of CDFTSVM  with ' num2str(frtsvm_struct.Parameter.ker),' kernel' ];
% title(line,'FontSize',12);

figure
sn=mapminmax(frtsvm_struct.sn',0.1,1)';
[xq,yq] = meshgrid(linspace(min(xn(:, 1)), max(xn(:, 1)), 100), linspace(min(xn(:, 2)), max(xn(:, 2)), 100));
vq = griddata(xn(:, 1), xn(:, 2),sn,xq,yq,'v4');
% vq(find(vq<0))=0;
contourf(xq,yq,vq,100,'LineStyle','none');
colormap(jet);
hold on
plot(xn(:, 1), xn(:, 2),'bx','LineWidth',1.5);
colorbar
hold off
% line=['Fuzzy membership distribution of CDFTSVM  with ' num2str(frtsvm_struct.Parameter.ker),' kernel' ];
% title(line,'FontSize',12);


figure
h1 = plot(xp(:, 1), xp(:, 2), 'r+','LineWidth',1.5);
text(xp(:,1),xp(:,2),num2str(round(100*sp)/100),'color','b','FontSize',4)%��(x,y)��дstring
hold on
h2 = plot(xn(:, 1), xn(:, 2), 'bx','LineWidth',1.5);
text(xn(:,1),xn(:,2),num2str(round(100*sn)/100),'color','r','FontSize',4)%��(x,y)��дstring
if ~isempty(frtsvm_struct.NXpv)|~isempty(frtsvm_struct.NXnv)
    if ~isempty(frtsvm_struct.NXpv)
    h3 = plot(xp(frtsvm_struct.NXpv,1),xp(frtsvm_struct.NXpv,2),'gs','MarkerSize',7);
    end
    if ~isempty(frtsvm_struct.NXnv)
    h3 = plot(xn(frtsvm_struct.NXnv,1),xn(frtsvm_struct.NXnv,2),'gs','MarkerSize',7);
    end
    if ~isempty(frtsvm_struct.NXpv)&&~isempty(frtsvm_struct.NXnv)
    h3 = plot(xp(frtsvm_struct.NXpv,1),xp(frtsvm_struct.NXpv,2),'gs','MarkerSize',7);    
    h3 = plot(xn(frtsvm_struct.NXnv,1),xn(frtsvm_struct.NXnv,2),'gs','MarkerSize',7);
    end
end
h4=legend([h1,h2,h3],'class +','class -','Outliers');
set(h4,'EdgeColor','w');
% line=['Fuzzy membership value of CDFTSVM  with ' num2str(frtsvm_struct.Parameter.ker),' kernel' ];
% title(line,'FontSize',12);
%%
figure;
clf;
set(gca,'XLim',[minX maxX],'YLim',[minY maxY]);
hold on;

h1 = plot(xp(:, 1), xp(:, 2), 'r+','LineWidth',1.5);
h2 = plot(xn(:, 1), xn(:, 2), 'bx','LineWidth',1.5 );
h3 = plot(psv(:,1),psv(:,2),'ko','MarkerSize',7 );
h3 = plot(nsv(:,1),nsv(:,2),'ko','MarkerSize',7 );
if ~isempty(frtsvm_struct.NXpv)|~isempty(frtsvm_struct.NXnv)
    if ~isempty(frtsvm_struct.NXpv)
    h4 = plot(xp(frtsvm_struct.NXpv,1),xp(frtsvm_struct.NXpv,2),'gs','MarkerSize',7);
    end
    if ~isempty(frtsvm_struct.NXnv)
    h4 = plot(xn(frtsvm_struct.NXnv,1),xn(frtsvm_struct.NXnv,2),'gs','MarkerSize',7);
    end
    if ~isempty(frtsvm_struct.NXpv)&&~isempty(frtsvm_struct.NXnv)
    h4 = plot(xp(frtsvm_struct.NXpv,1),xp(frtsvm_struct.NXpv,2),'gs','MarkerSize',7);    
    h4 = plot(xn(frtsvm_struct.NXnv,1),xn(frtsvm_struct.NXnv,2),'gs','MarkerSize',7);
    end
end

h4=legend([h1,h2,h3,h4],'class +','class -','Support Vectors','Outliers');
set(h4,'EdgeColor','w');
[C,h] = contour(bigX, bigY, bigZ,[0 0],'k','LineWidth',1.5);
clabel(C,h,'Color','k');
[C,h] = contour(bigX, bigY, bigf1,[0 0],'b:','LineWidth',1.5);
clabel(C,h,'Color','b');
[C,h] = contour(bigX, bigY, bigf2,[0 0],'r:','LineWidth',1.5);
clabel(C,h,'Color','r');

xlabel('demension1','FontSize',12);
ylabel('demension2','FontSize',12);
% line=['Visualization of',' CDFTSVM classification with ' num2str(frtsvm_struct.Parameter.ker),' kernel' ];
% title(line,'FontSize',12);
end

