% f=-[3;4;7;8];
% A=[10,15,25,30];
% b=50;
% lb=zeros(n,1);
% ub=ones(n,1);
% [x,fval,exitflag,OUTPUT] = intlinprog(f,1:4,A,b,[],[],lb,ub,[]);
for s=2:62
    sprintf("s:%d, fix0:%d, fix1:%d, nearzero_intgap_itr:%d", s, length(fix63(s).fixto0list), ...
        length(fix63(s).fixto1list), fix63(s).nearzero_intgap_itr)
    if length(fix63(s).fixto0list)+length(fix63(s).fixto1list) == 63
        s
    end
end