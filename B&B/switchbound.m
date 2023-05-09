% Switch type of bound (original/complement) and check if other bound
% gives better value.  If so update bound and switch value of "complement"
%
complement = ~complement;                       % switch
control(6)=0;                                   % re-initialize gamma value
if complement                                   % apply bound to complement of node problem
    [results2,xval2,delta_one2,delta_zero2]=eval(strcat(solver,'(Cnodeinv,nnode-snode,control)')); % note deltas switched            
    nodeconstant2=fixconstant+logdetCnode;
    xval2=ones(nnode,1)-xval2;
else
    [results2,xval2,delta_zero2,delta_one2]=eval(strcat(solver,'(Cnode,snode,control)'));           % use original problem
    nodeconstant2=fixconstant;
end     
bound2=results(1)+nodeconstant2;
code=results(2);
solvetime=solvetime+results(3);
kscale=results(5);
Ncompute2=Ncompute2+1;
%
Nscale(kscale)=Nscale(kscale)+1;
%
if code < 0
    Ncodeneg=Ncodeneg+1;
elseif code == 0
    Ncode0=Ncode0+1;
else Ncodepos(code)=Ncodepos(code)+1;
end
%
if bound2 < fathomval       % node is fathomed
    Nbetter2=Nbetter2+1;    % second bound was better
    clearnode;              % done with node
    continue 
elseif bound2 < bound;      % switch to results from better bound
    bound=bound2;
    Nbetter2=Nbetter2+1;    % second bound was better
    nodeconstant=nodeconstant2;
    delta_zero=delta_zero2;
    delta_one=delta_one2;
    xval=xval2;
    newgamma=results(4);
else complement = ~complement; % go back - no reason to switch
end