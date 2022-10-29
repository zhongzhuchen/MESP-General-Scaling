% B&B driver for MESP solvers.  Considers dynamic use of 
% original and complimentary problems. Uses initial dive to attempt
% to find optimal solution earlier. Checkpoints results to allow restart.
% Tolerances for fathoming and fixing variables can be increased to 
% positive values after first checkpoint.
%
clear all;
%
restart = 0;    % 1 if restarting from checkpoint
tic; % start timer
%                             
if restart      % restart from checkpoint file
    load('checkpoint1.mat');
else         
    %      ******* problem initialization ********
    load('data124.mat');
    n=length(C);
%     fid = fopen('DATA63','r');
%     %
%     [C,count] = fscanf(fid,'%f',[n,n]);
%     fclose(fid);
    %
    s=30;
    complement=0;       % initial setting for bound (original or complement)
    %
    [U,D]=eig(C);
    lam=diag(D);
    logdetC=log(prod(lam));  
    Cinv=U*diag(1./lam)*U';
    %
    [xind,heurval]=heur1(C,n,s);        % HEURSITIC ON ORIGINAL
    [cxind,cheurval]=heur1(Cinv,n,n-s); % HEURISTIC ON COMPLEMENT
    if cheurval+logdetC > heurval      % PICK THE BEST
      xind=setdiff(transpose(1:n),cxind(1:n-s));
      heurval=cheurval+logdetC;
    end;
    fprintf('\n HEURISTIC VALUE=%g\n\n', heurval)
    %
    %  ***** Tolerances for fathoming nodes and fixing variables  ******
    %
    fathomtol=0;      % tolerance for fathoming. If > 0 fathom nodes within tol of bestval
    fathomtol2=0;     % change to fathomtol2 after first checkpoint
    fixtol=0;         % tolerance for fixing variables. If fixtol > 0 then
                      % variable-fixing logic will eliminate solutions within
                      % fixtol of bestval. Used in fixvariables 
    fixtol2=0;        % change to fixtol2 after first checkpoint      
    besttol=-.1;      % tolerance for initial best known value. 
                      % set < 0 to insure generation of optimal solution
    %    
    bestval=heurval+besttol;
    fathomval=bestval+fathomtol; % adjust fathoming criterion using tolerance
    fixval=bestval+fixtol;
    bestvars=zeros(1,s);   % indeces of variables for bestval

    Qprob=zeros(500,6-1+n);    % problem queue: Parent Bound,nfix0,nfix1,depth,orig/comp,newgamma
    % Todo, the problem queue may overflow
    Qprob(1,1)=10000;      % prevent fathoming based on parent bound
    Qprob(1,5)=complement; % type of bound
    Qfix0=zeros(500,n-s);  % variables fixed to zero in problem queue
    Qfix1=zeros(500,s);    % variables fixed to one in problem queue
                           % assume indeces for fixed variables are in first
                           % nfix0/nfix1 positions
    %
    nq=1;           % number of nodes in queue
    Nnodes=0;       % number of nodes processed
    Nfathomed=0;    % number of nodes fathomed    
    Nenumerate=0;   % number fathomed by enumeration  
    Nfractional=0;  % number of branchings based on fractional variables
    Ncompute2=0;    % number of times second bound computed
    Nbetter2=0;     % number of times second bound better
    Ndepth=zeros(n,1);  % Ndepth(i+1)= number of nodes at depth i in tree (root is depth 0)
    Nfathom=zeros(n,1); % number of nodes fathomed at depth i
    Nnfix=zeros(n,1);   % keep track of sizes of subproblems
    Ncodeneg=0;         % number of negative solver codes
    Ncode0=0;           % number of zero solver codes
    Ncodepos=zeros(10,1); % vector for positive solver codes
    Nscale=zeros(10,1); % number of scalings used in bound comps
    Nscalerr0=0;        % save results of scale errors from bound comps (log2): error in (.5, 1]
    Nscalerrneg=zeros(10,1);    % -1: error in (.25,.5]. -2: error in [.125, .25]. etc.
    Nscalerrpos=zeros(10,1);
    Qmax=1;             % max length of queue
    ndive=1;            % last node in tree from initial dive
    time=0;             % total run time
    solvetime=0;        % time spent in solvers
    toggle=0;           % to alternate checkpoint files
end
%
% ***** B&B Settings - note can change settings after restart ******
%
freq=5;             % consider switching orig/comp if mod(depth,freq)=0 (never if > n)
frac=1;             % and at least frac of required # of 0 or 1 variables are fixed (never if =1)
deltaminabs=.025;   % absolute tolerance for branching using dual info
deltaminfrac=.1;    % relative tolerance factor for branching using dual info
minfractional=0;    % switch to fractional branching if nfractional <= minfractional (0 for never)
fractol=.005;       % tolerance to consider variable to be fractional 
%
% **** Solver and control variables for solver *******
%
solver='linx';
%
nscale=2;     % max number of scale facors to usually try
nscale0=5;    % max number on initial bound evaluation
nrounds=0;    % number of rounds adding inequalities (not used by linx)
Maxadd=2000;  % max number of inequalities per round (not used by linx)
printsol=0;   % print solver results
Gamma=0;      % initial gamma value
power=1.5;    % power for initial gamma computation
maxfactor=5;  % max factor change for gamma on one update
%                  
% control=[nscale, nscale0, nrounds, Maxadd, printsol, gamma, power, maxfactor];
control=struct;
control.nscale = nscale;
control.nscale0 = nscale0;
control.nrounds = nrounds;
control.Maxadd = Maxadd;
control.printsol = printsol;
control.Gamma = Gamma;
control.power = power;
control.maxfactor = maxfactor;

while nq > 0              % queue is not empty
    elapsedtime=toc;
    if elapsedtime > 1800 % interval for checkpointing
        time=time+elapsedtime;
        if toggle; save('checkpoint1.mat');
        else save('checkpoint0.mat');
        end
        fathomtol=fathomtol2;  % update tolerances
        fixtol=fixtol2;
        fathomval=bestval+fathomtol;
        fixval=bestval+fixtol;
        %
        toggle = 1-toggle; % alternate checkpoint file
        tic; % restart timer
    end
    %
    % Take last node off queue
    %
    Nnodes=Nnodes+1;
    if nq > Qmax; Qmax = nq; end
    bound=Qprob(nq,1);      % bound from parent problem
    nfix0=Qprob(nq,2);      % number of variables fixed to zero
    nfix1=Qprob(nq,3);      % number of variables fixed to one
    depth=Qprob(nq,4);      % depth in B&B tree
    complement=Qprob(nq,5); % original or complimentary problem for parent bound
    %
    Ndepth(depth+1) = Ndepth(depth+1)+1; % note root is at depth 0
    % 
    if nq==ndive            % node is on inital dive going down, or child of dive node coming back up
        fprintf('nq=%g, Nnodes=%g, Bound=%g, Bestval=%g, Time=%g, nfix0=%g,nfix1=%g\n',nq,Nnodes,bound,bestval,time+toc,nfix0,nfix1); 
        if Nfathomed == 0   % going down
            ndive=ndive+1;
        else ndive=ndive-1; % coming back up
        end
    end
    %
    vfix0=Qfix0(nq,:); % note fixed variables are in row vectors
    vfix1=Qfix1(nq,:);
    %
    % Error logic to ttop if repeated values; note also count "0"
    % if size(unique(vfix0),2) ~= nfix0+1 || size(unique(vfix1),2) ~= nfix1+1 break; end
    %   
    vfix=union(vfix0,vfix1);        % includes "0" - cannot have all fixed
    vars=setdiff(0:n,vfix);         % remaining variables for node subproblem
    nfix=nfix0+nfix1;
    nnode=n-nfix; 
    snode=s-nfix1;
    Nnfix(nfix+1)=Nnfix(nfix+1)+1;  % start with nfix=0 at root
    %
    if bound < fathomval            % check parent bound; fathomval may have changed
        clearnode;                  % done with node
        continue                    
    end    
    %
    if snode <=1 || nnode-snode <=1
        enumerate;  % enumerate solutions
        clearnode;  % done with node
        continue 
    end
    %   
    Cnode=C(vars,vars); %submatrix for node
    fixconstant=0;    
    if nfix1>0
        vvfix1=vfix1(1:nfix1);      % cut vector down to actual fixed components
        [U,D]=eig(C(vvfix1,vvfix1));
        lam=diag(D);
        fixconstant=log(prod(lam)); % adjustment for fixed variables 
        Cfixinv=U*diag(1./lam)*U';
        Cnode=Cnode-C(vars,vvfix1)*Cfixinv*C(vvfix1,vars); %Shur complement
        Cnode=.5*(Cnode+Cnode');    % force symmetrization to avoid errors
    end
    [U,D]=eig(Cnode);
    lam=diag(D);
    logdetCnode=log(prod(lam));  
    Cnodeinv=U*diag(1./lam)*U';     % Compute inverse for use in complementary bound
    % take Gamma according to vars
    Gammaind = vars+5;
    Gamma=Qprob(nq,Gammaind);      % final adjusted gamma from parent problem
    control.Gamma=Gamma;               % use final adjusted gamma from parent problem
    if complement                   % apply bound to complement of node problem
        [results,xval,delta_one,delta_zero]=eval(strcat(solver,'(Cnodeinv,nnode-snode,control)')); % note deltas switched  
        nodeconstant=fixconstant+logdetCnode;
        xval=ones(nnode,1)-xval;
    else 
        [results,xval,delta_zero,delta_one]=eval(strcat(solver,'(Cnode,snode,control)'));           % use original problem
        nodeconstant=fixconstant;
    end
    %
    bound=results.bound+nodeconstant;  % record outputs from bound computation
    code=results.code;
    solvetime=solvetime+results.time;
    newGamma=results.newGamma;
    newGammafull = ones(1,n);
    newGammafull(vars) = newGamma; % supplement newGamma to of length n for inheritance
    kscale=results.kscale;
    scalerror=results.scalerror;
    %
    if code < 0                     % solver solution code
        Ncodeneg=Ncodeneg+1;
    elseif code == 0
        Ncode0=Ncode0+1;
    else Ncodepos(code)=Ncodepos(code)+1;
    end
    %
    Nscale(kscale)=Nscale(kscale)+1;    % number of gamma scale factors used
    lscalerror=ceil(log2(scalerror));   % record error in scaling criterion (log2)
    lscalerror=min(lscalerror,10);
    lscalerror=max(lscalerror,-10);
    if lscalerror > 0 
        Nscalerrpos(lscalerror)=Nscalerrpos(lscalerror)+1;
    elseif lscalerror < 0
        Nscalerrneg(abs(lscalerror))=Nscalerrneg(abs(lscalerror))+1;
    else Nscalerr0=Nscalerr0+1;
    end             
    %
    if bound < fathomval    % node is fathomed
        clearnode;          % done with node
        continue 
    end
    %
    checkfrac = nfix1 > frac*s || nfix0 > frac*(n-s);   % criteria for number of fixed variables
    if mod(depth,freq)==0 && checkfrac                  % consider swithing orig/comp problems
        swithbound;
    end
    %
    fixvariables;   % use dual info to fix variables to 0/1 values
    %
    if snode <=1 || nnode-snode <=1
        enumerate;  % enumerate solutions
        clearnode;  % done with node
        continue 
    end
    %
    nfractional=0;
    for i=1:nnode;  % check number of fractional variables
        if xval(i) > fractol && xval(i) < 1-fractol;
            nfractional=nfractional+1;
        end
    end
    %
    fancybranch;    % generate child problems
end
%
time=time+toc;
fprintf('\n Nnodes=%g, Overall time=%g, Solver time=%g, Bestval=%g\n', Nnodes, time, solvetime, bestval);
    
