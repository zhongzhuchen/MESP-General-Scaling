% Clear fathomed node from queue
%
Qprob(nq,:)=zeros(1,6);
Qfix0(nq,:)=zeros(1,n-s);
Qfix1(nq,:)=zeros(1,s);
%
Nfathom(depth+1)=Nfathom(depth+1)+1; % record node as fathomed
if Nfathomed == 0
    ndive=nq-1; % set end of initial dive
end
Nfathomed=Nfathomed+1;
nq=nq-1;
