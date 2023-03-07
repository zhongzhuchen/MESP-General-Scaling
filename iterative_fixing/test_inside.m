load('fix300_80to120_LFcF_iter1.mat');
fixnumsnew1 = fixnumsnew;
load('fix300_80to120_LFcF_iter2.mat');
fixnumsnew2 = fixnumsnew;
load('fix300_80to120_LFcF_iter3.mat');
fixnumsnew3 = fixnumsnew;
load('fix300_80to120_LFcF_iter4.mat');
fixnumsnew4 = fixnumsnew;
load('fix300_80to120_LFcF_iter5.mat');
fixnumsnew5 = fixnumsnew;
load('fix300_80to120_LFcF_iter10.mat');
fixnumsnew10 = fixnumsnew;
load('fix300_80to120_LFcF_iter20.mat');
fixnumsnew20 = fixnumsnew;
load('fix300_80to120_LFcF_iter100.mat');
fixnumsnew100 = fixnumsnew;

fixnums_compare = [fixnumsnew1.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew2.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew3.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew4.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew5.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew10.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew20.gscaling_300'-fixnumsnew1.oscaling_300',...
    fixnumsnew100.gscaling_300'-fixnumsnew1.oscaling_300']
% fixnumsnew = struct;
% clearvars -except fixnumsnew fixnums;
% load('data63.mat');
% iterative_fixing_G;
% iterative_fixing_O;
% clearvars -except fixnumsnew fixnums;
% load('data90.mat');
% iterative_fixing_G;
% iterative_fixing_O;
% clearvars -except fixnumsnew fixnums;
% load('data124.mat');
% iterative_fixing_G;
% iterative_fixing_O;
 
%==========================================================================
% load('fix300_80to120_LFcF_gscaling.mat');
% load('data2000_sub300.mat');
% n=length(C);
% s=105;
% vars = setdiff(1:n, union(fix(s).fixto0list, fix(s).fixto1list));
% vars = sort(vars);
% Cnode=C(vars,vars);
% vvfix1=fix(s).fixto1list;
% [U,D]=eig(C(vvfix1,vvfix1));
% lam=diag(D);
% fixconstant=log(prod(lam)); % adjustment for fixed variables 
% Cfixinv=U*diag(1./lam)*U';
% Cnode=Cnode-C(vars,vvfix1)*Cfixinv*C(vvfix1,vars); %Shur complement
% Csub=.5*(Cnode+Cnode'); 
%==========================================================================
% save('Csub.mat', "Csub");
% load('data63.mat');
% constr_iterative_fixing_G;
% load('data90.mat');
% constr_iterative_fixing_G;
% load('data124.mat');
% constr_iterative_fixing_G;
%==========================================================================
% solved_instance = struct;
% for n=[63]
%     filename = strcat('constr_fix', num2str(n),'_LFcF_oscaling.mat');
%     load(filename);
%     % load('fix63_LFcF_oscaling.mat');
%     solvelist = [];
%     fixalllist1 = [];
%     fixnum1 = [];
%     for s=3:52
%         if fix(s).nearzero_intgap == 1
%             solvelist(end+1)=s;
%         end
%         if length(fix(s).fixto0list)+length(fix(s).fixto1list) >= n
%             fixalllist1(end+1)=s;
%         end
%         fixnum1(end+1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
%     end
% %     fixalllist1 = union(solvelist, fixalllist1)
% %     solved_instance.oscaling_63_solvelist = fixalllist1;
%     if n==63
%         fixnums.oscaling_63 = fixnum1;
%     elseif n==90
%         fixnums.oscaling_90 = fixnum1;
%     elseif n==124
%         fixnums.oscaling_124 = fixnum1;
%     end
% 
%     filename = strcat('constr_fix', num2str(n),'_LFcF_gscaling.mat');
%     load(filename);
%     % load('fix63_LFcF_gscaling.mat');
%     solvelist = [];
%     fixalllist2 = [];
%     fixnum = [];
%     for s=3:52
%         if fix(s).nearzero_intgap == 1
%             solvelist(end+1)=s;
%         end
%         if length(fix(s).fixto0list)+length(fix(s).fixto1list) >= n
%             fixalllist2(end+1)=s;
%         end
%         fixnum(end+1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
%     end
% %     fixalllist2 = union(solvelist, fixalllist2)
% %     solved_instance.gscaling_63_solvelist = fixalllist2;
%     if n==63
%         fixnums.gscaling_63 = fixnum;
%     elseif n==90
%         fixnums.gscaling_90 = fixnum;
%     elseif n==124
%         fixnums.gscaling_124 = fixnum;
%     end
% %     [c,~] = ismember(fixalllist1,fixalllist2)
% end

%==========================validality check================================
% load('data63.mat');
% load('linconstr63.mat');
% load('constr_fix63_LFcF_gscaling.mat');
% MESPInstance = MESP(C,A,b);
% F=MESPInstance.F;
% Fsquare=MESPInstance.Fsquare;
% F_comp=MESPInstance.F_comp;
% Fsquare_comp=MESPInstance.Fsquare_comp;
% ldetC=MESPInstance.ldetC;
% n=length(C);
% Gamma=ones(n,1);
% ubcode=[];
% fixnums = [];
% for s=3:52
%     x0=s/n*ones(n,1);
%     [fval,x,info1] = Knitro_Linx_light_check(x0,C,s,A,b,Gamma,fix(s));
%     [fval,x,info2]=Knitro_DDFact_light_check(x0,C,s,F,Fsquare,A,b,Gamma,fix);
%     [fval,x,info3] = Knitro_DDFact_comp_light_check(x0,C,s,F_comp,Fsquare_comp,ldetC,A,b,Gamma,fix);
%     [xind, heurval] = heur(C,s,A,b);
%     ubcode(end+1,:)=[info1.fval-heurval>=0,info2.fval-heurval>=0,info3.fval-heurval>=0,info.exitflag];
%     
% %     [xind] = gencontsol_diag_check(C,s,A,b,fix);
% end







