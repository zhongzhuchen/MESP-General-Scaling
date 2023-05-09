% Branching logic for MESPBB
% uses info from fixvariables
%
branchfractional=1;              % logic to use most fractional variable versus dual information
mindelta=min(mindelta_zero,mindelta_one);
if  nfractional > minfractional  % otherwise branch on fractional variable regardless of deltas
    if mindelta < -deltaminabs || mindelta < deltaminfrac*(fixval-bound) % criteria for branching based on mindelta
        branchfractional=0;
    end
end
%
if Nfathomed == 0; % still on initial dive; put harder branch on last
   if mindelta_zero < mindelta_one % will branch on branchvar_zero
        % start with easier node: branchvar_zero=0
        Qprob(nq,1)=bound+mindelta_zero;
        Qprob(nq,2)=nfix0+1; 
        Qprob(nq,3)=nfix1;
        Qprob(nq,4)=depth+1;
        Qprob(nq,5)=complement;
        Qprob(nq,6)=newgamma;
        %
        vfix0(nfix0+1)=branchvar_zero;
        Qfix0(nq,:)=vfix0;
        Qfix1(nq,:)=vfix1;
        vfix0(nfix0+1)=0; % go back
        % now add node with branchvar_zero=1
        nq=nq+1;
        Qprob(nq,1)=bound;
        Qprob(nq,2)=nfix0;
        Qprob(nq,3)=nfix1+1;
        Qprob(nq,4)=depth+1;
        Qprob(nq,5)=complement;
        Qprob(nq,6)=newgamma;
        %
        vfix1(nfix1+1)=branchvar_zero;
        Qfix1(nq,:)=vfix1;
        Qfix0(nq,:)=vfix0;
        vfix1(nfix1+1)=0; % go back
   else % will branch on branchvar_one 
        % start with easier node: branchvar_one=1
        Qprob(nq,1)=bound+mindelta_one;
        Qprob(nq,2)=nfix0;
        Qprob(nq,3)=nfix1+1;
        Qprob(nq,4)=depth+1;
        Qprob(nq,5)=complement;
        Qprob(nq,6)=newgamma;
        %
        vfix1(nfix1+1)=branchvar_one;
        Qfix1(nq,:)=vfix1;
        Qfix0(nq,:)=vfix0;
        vfix1(nfix1+1)=0; % go back
        % now add node with branchvar_one=0
        nq=nq+1;
        Qprob(nq,1)=bound;
        Qprob(nq,2)=nfix0+1; 
        Qprob(nq,3)=nfix1;
        Qprob(nq,4)=depth+1;
        Qprob(nq,5)=complement;
        Qprob(nq,6)=newgamma;
        %
        vfix0(nfix0+1)=branchvar_one;
        Qfix0(nq,:)=vfix0;
        Qfix1(nq,:)=vfix1;
        vfix0(nfix0+1)=0; % go back
   end
elseif branchfractional % branch on most fractional variable; order does not matter
            Nfractional=Nfractional+1;
            Qprob(nq,1)=bound;
            Qprob(nq,2)=nfix0;
            Qprob(nq,3)=nfix1+1;
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix1(nfix1+1)=maxfracvar;
            Qfix1(nq,:)=vfix1;
            Qfix0(nq,:)=vfix0;
            vfix1(nfix1+1)=0; % go back
            % 
            nq=nq+1;
            Qprob(nq,1)=bound;
            Qprob(nq,2)=nfix0+1; 
            Qprob(nq,3)=nfix1;
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix0(nfix0+1)=maxfracvar;
            Qfix0(nq,:)=vfix0;
            Qfix1(nq,:)=vfix1;
            vfix0(nfix0+1)=0; % go back    
elseif nhalf <= 1 % Dive is over, but not enough half variables to do anything special
                  % Branch as before, but switch order to put harder node on first
       if mindelta_zero < mindelta_one % will branch on branchvar_zero
            % start with branchvar_zero=1
            Qprob(nq,1)=bound;
            Qprob(nq,2)=nfix0;
            Qprob(nq,3)=nfix1+1;
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix1(nfix1+1)=branchvar_zero;
            Qfix1(nq,:)=vfix1;
            Qfix0(nq,:)=vfix0;
            vfix1(nfix1+1)=0; % go back
            % now add easier node: branchvar_zero=0
            nq=nq+1;
            Qprob(nq,1)=bound+mindelta_zero;
            Qprob(nq,2)=nfix0+1; 
            Qprob(nq,3)=nfix1;
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix0(nfix0+1)=branchvar_zero;
            Qfix0(nq,:)=vfix0;
            Qfix1(nq,:)=vfix1;
            vfix0(nfix0+1)=0; % go back
       else % will branch on branchvar_one; 
            % start with branchvar_one=0
            Qprob(nq,1)=bound;
            Qprob(nq,2)=nfix0+1; 
            Qprob(nq,3)=nfix1;
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix0(nfix0+1)=branchvar_one;
            Qfix0(nq,:)=vfix0;
            Qfix1(nq,:)=vfix1;
            vfix0(nfix0+1)=0; % go back
            % now add node with branchvar_one=1
            nq=nq+1;
            Qprob(nq,1)=bound+mindelta_one;
            Qprob(nq,2)=nfix0;
            Qprob(nq,3)=nfix1+1;
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix1(nfix1+1)=branchvar_one;
            Qfix1(nq,:)=vfix1;
            Qfix0(nq,:)=vfix0;
            vfix1(nfix1+1)=0; % go back           
       end
elseif nfix1+nhalf0 > s || nfix0+nhalf1 > n-s % have enough half variables to get fancy
       clearnode;                             % first check if branching is possible
else vfix0old=vfix0; % save for re-use
     vfix1old=vfix1;
    %
    %first put on node that doesn't use any half variables
        Qprob(nq,1)=bound;
        Qprob(nq,2)=nfix0+nhalf1; % all half1 variables fixed at zero
        Qprob(nq,3)=nfix1+nhalf0; % all half0 variables fixed at 1
        Qprob(nq,4)=depth+1;
        Qprob(nq,5)=complement;
        Qprob(nq,6)=newgamma;
        %
        if nhalf1 > 0
            vfix0(nfix0+1:nfix0+nhalf1)=vhalf1(1:nhalf1);
        end
        %
        if nhalf0 > 0
            vfix1(nfix1+1:nfix1+nhalf0)=vhalf0(1:nhalf0);
        end
        Qfix1(nq,:)=vfix1;
        Qfix0(nq,:)=vfix0;
    %
    if nfix0+nhalf1+1 <= n-s % check if possible first
        for i=1:nhalf0 % add one node for each of the half0 variables
            nq=nq+1;   
            vfix0=vfix0old;
            vfix1=vfix1old;
            %
            Qprob(nq,1)=bound; % note don't currently have delta info for bound
            Qprob(nq,2)=nfix0+nhalf1+1; % all half1 variables fixed at zero
            Qprob(nq,3)=nfix1+nhalf0-1; % all other half0 variables fixed at 1
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix0(nfix0+1)=vhalf0(i);
            if nhalf1 > 0
                vfix0(nfix0+2:nfix0+nhalf1+1)=vhalf1(1:nhalf1);
            end
            if nhalf0 > 1
               others=setdiff((1:nhalf0),i);
               vfix1(nfix1+1:nfix1+nhalf0-1)=vhalf0(others);
            end 
            Qfix1(nq,:)=vfix1;
            Qfix0(nq,:)=vfix0;
        end
   end
    %
    if nfix1+nhalf0+1 <= s % check if possible first
        for i=1:nhalf1 % add one node for each of the half1 variables
            nq=nq+1;
            vfix0=vfix0old;
            vfix1=vfix1old;
            %
            Qprob(nq,1)=bound;
            Qprob(nq,2)=nfix0+nhalf1-1; % all other half1 variables fixed at zero
            Qprob(nq,3)=nfix1+nhalf0+1; % all half0 variables fixed at 1
            Qprob(nq,4)=depth+1;
            Qprob(nq,5)=complement;
            Qprob(nq,6)=newgamma;
            %
            vfix1(nfix1+1)=vhalf1(i);
            if nhalf0 > 0
                vfix1(nfix1+2:nfix1+nhalf0+1)=vhalf0(1:nhalf0);
            end
            if nhalf1 > 1
               others=setdiff((1:nhalf1),i);
               vfix0(nfix0+1:nfix0+nhalf1-1)=vhalf1(others);
            end 
            Qfix1(nq,:)=vfix1;
            Qfix0(nq,:)=vfix0;
        end
    end 
    vfix0=vfix0old;
    vfix1=vfix1old;
end


       
       
       
       
      
