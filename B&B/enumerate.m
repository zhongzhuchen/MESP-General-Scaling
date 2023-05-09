% enumerate solutions for MESPBB if snodes<=1 or nnodes-snodes<=1
% Note uses bestval as criterion.
%
oldbest=bestval;
Nenumerate=Nenumerate+1;
%
if snode==0 % nfix1=s => no more left to choose
    value=log(det(C(vfix1,vfix1)));
    if value>bestval;
        bestval=value;
        bestvars=vfix1;       
    end
elseif snode==nnode % nfix0 = n-s => must use all the rest
    vfix1=setdiff(1:n,vfix0);
    value=log(det(C(vfix1,vfix1)));
    if value>bestval;
        bestval=value;
        bestvars=vfix1;
    end
elseif snode==1; % s - nfix1 = 1 => must add one more variable
    for i=1:nnode
        vfix1(s)=vars(i); % choose last variable
        value=log(det(C(vfix1,vfix1)));
        if value>bestval;
            bestval=value;
            bestvars=vfix1;
        end
    end
elseif nnode-snode == 1 % n-s = nfix0 + 1 => must set one more to zero
    for i=1:nnode;
        vfix0(n-s)=vars(i); %choose last to set to zero
        vfix1=setdiff(1:n,vfix0);
        value=log(det(C(vfix1,vfix1)));
        if value>bestval;
            bestval=value;
            bestvars=vfix1;   
        end
    end
end
%
if bestval > oldbest; % found new best value
    fathomval=bestval+fathomtol;
    fixval=bestval+fixtol;
end
