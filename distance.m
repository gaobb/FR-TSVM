function [di,xstar]=distance(xi,X)
%  Author: Bin-Bin Gao
%  Email: csgaobb@gmail.com
%  July 5, 2017

% check correct number of arguments
if ( nargin>2||nargin<2)
    help distance
else
    l =size(X,1);
    for  i=1:l
        di(i,1)=norm(xi-X(i,:));
    end
    sx =X(min(di) ==di,:);
    xstar = sx(1,:);
end
end
