% variable-fixing logic for MESPBB
% Note fixing logic is based on fixval.
%
mindelta_zero=1; % min of deltas not good enough to fix variable to 0/1
mindelta_one=1;
nhalf0=0;
nhalf1=0;
vhalf0=zeros(1,nnode);
vhalf1=zeros(1,nnode);
maxfraction=0; % largest fractional part
%
for i=1:nnode
    if bound+delta_one(i)< fixval % can fix variable to zero
        nfix0=nfix0+1;
        vfix0(nfix0)=vars(i);
    elseif bound+delta_zero(i)< fixval % can fix variable to one
        nfix1=nfix1+1;
        vfix1(nfix1)=vars(i);        
    else % variable cannot be fixed
        if delta_one(i) < (fixval-bound)/2 && delta_one(i) >= fixval-bound
            nhalf1=nhalf1+1;  % variables with delta_one halfway to fixing value 
            vhalf1(nhalf1)=vars(i); 
        elseif delta_zero(i) < (fixval-bound)/2 && delta_zero(i) >= fixval-bound
            nhalf0=nhalf0+1; % variables with delta_zero halfway to fixing value
            vhalf0(nhalf0)=vars(i);
        end
        %
        if delta_one(i) < mindelta_one % check for best remaining
            mindelta_one=delta_one(i);
            branchvar_one=vars(i);
        end
        if delta_zero(i) < mindelta_zero % check for best remaining
            mindelta_zero=delta_zero(i);
            branchvar_zero=vars(i);
        end
        %
        if min(xval(i),1-xval(i)) > maxfraction
            maxfraction = min(xval(i), 1-xval(i));
            maxfracvar=vars(i);
        end
    end
end
%
vfix=union(vfix0,vfix1);
nfix=nfix0+nfix1;
nhalf=nhalf0+nhalf1;
if nfix<n; vars=setdiff(0:n,vfix); end % remaining variables; allow for all fixed
%
nnode=n-nfix;
snode=s-nfix1;
%








            
            
        
    