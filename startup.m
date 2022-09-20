load('data63.mat');
n=length(C);
% load('linconst63.mat')
A2 = double.empty(0,n);
b2 = double.empty(0,1);

% for s=5:n-15
%     [xind, heurval] = heur(C,s,A2,b2);
%     x=zeros(n,1);
%     x(xind)=1;
%     xlist(:,s-4)=A*x;
% end
% % sort by row
% xlist=sort(xlist,2);
% ind=round((n-9)/5*4);
% b=xlist(:,ind)-1;
% anchor=1:8:39;
% anchor(6)=40;
% b=zeros(5,1);
% for i=1:5
%     range=anchor(i):(anchor(i+1)-1);
%     b(i)=min(xlist(i,range))-1;
% end
MESPInstance = MESP(C,A,b);

MESPInstance2 = MESP(C,A2,b2);