clear all;
load('data63.mat');
n=length(C);
folder = 'old_Results';
baseFileName = strcat('data',int2str(n),'_mix_uniform_lincon','.xlsx');
fullFileNameexcel = fullfile(folder, baseFileName);
if exist(fullFileNameexcel, 'file')==2
  delete(fullFileNameexcel);
end

A2 = double.empty(0,n);
b2 = double.empty(0,1);

A = randi([-2,2],5,n);

xlist=[];
% generate linear constraint
deduct = 13;
start = 5;
for s=start:(n-deduct)
    [xind, heurval] = heur(C,s,A2,b2);
    x=zeros(n,1);
    x(xind)=1;
    xlist(:,s-start+1)=A*x;
end
% xlist=xlist(:, randperm(size(xlist, 2)));
% sort by row
% xlist=sort(xlist,2);
% ind=round((n-9)/5*4);
% b=xlist(:,ind)-1;
anchor=1:8:(n-deduct-start+1);
b=zeros(5,1);
for i=1:5
    range=anchor(i):(anchor(i+1)-1);
    b(i)=min(xlist(i,range))-1;
end
% for i=1:5
%     b(i)=min(xlist(i,anchor(i)))-1;
% end
MESPInstance =  MESP(C,A,b);
MESPInstance2 = MESP(C,A2,b2);
% Prob = MESP(C,A,b);
% totaltimes=1;
% for repeattime = 1:totaltimes
%     exceloutput=[];
%     for s=5:(n-16)
%         sprintf("In the current iteration, n=%d, s=%d, repeattime=%d", n, s, repeattime)
%         %% initialize the MESP instance for each (m,s) pair
%         
%         lb = Prob.obtain_lb(s);
%     
%         x0=s/n*ones(n,1);
%         %% Linx data
%         [gamma,info]= Prob.BFGS_Linx_gamma(s);
%         [~,~,info_Linx] = Prob.Knitro_Linx(x0,s,sqrt(gamma)*ones(n,1));
%        
%         %% DDFact data
%         [~,~,info_DDFact] = Prob.Knitro_DDFact(x0,s,ones(n,1));
%     
%         %% DDFact comp data
%         [~,~,info_DDFact_comp] = Prob.Knitro_DDFact_comp(x0,s,ones(n,1));
%         
% %         info_mix1.dualbound=nan;
% %         info_mix1.integrality_gap=nan;
% %         info_mix1.continuous_dualgap=nan;
% %         info_mix1.fixnum_to0=nan;
% %         info_mix1.fixnum_to1=nan;
% %         info_mix1.fixnum=nan;
% %         info_mix1.alpha=nan;
% %         info_mix1.time=nan;
% %         info_mix1.cputime=nan;
% % 
% %         info_mix2.dualbound=nan;
% %         info_mix2.integrality_gap=nan;
% %         info_mix2.continuous_dualgap=nan;
% %         info_mix2.fixnum_to0=nan;
% %         info_mix2.fixnum_to1=nan;
% %         info_mix2.fixnum=nan;
% %         info_mix2.alpha=nan;
% %         info_mix2.time=nan;
% %         info_mix2.cputime=nan;
% % 
% %         info_mix3.dualbound=nan;
% %         info_mix3.integrality_gap=nan;
% %         info_mix3.continuous_dualgap=nan;
% %         info_mix3.fixnum_to0=nan;
% %         info_mix3.fixnum_to1=nan;
% %         info_mix3.fixnum=nan;
% %         info_mix3.alpha=nan;
% %         info_mix3.time=nan;
% %         info_mix3.cputime=nan;
%         if repeattime == 1
%             %% mix DDFact and/or comp
%             [~,~,info_mix1] = Prob.mix_DDFact_DDFact_comp(x0,s,ones(n,1),ones(n,1));
%             %% mix DDFact Linx
%             [optgamma, ~] = Prob.mix_BFGS_Gamma(s, "DDFact_Linx", gamma);
%             [~,~,info_mix2] = Prob.mix_DDFact_Linx(x0,s,ones(n,1),sqrt(optgamma)*ones(n,1));
%         
%             %% mix DDFact comp Linx
%             [optgamma, ~] = Prob.mix_BFGS_Gamma(s, "DDFact_comp_Linx", gamma);
%             [~,~,info_mix3] = Prob.mix_DDFact_comp_Linx(x0,s,ones(n,1),sqrt(optgamma)*ones(n,1));
%         end
% 
%         exceloutput(end+1,:)=[n, s, lb,...
%             info_Linx.exitflag, info_DDFact.exitflag, info_DDFact_comp.exitflag,...
%             info_Linx.dualbound, info_DDFact.dualbound, info_DDFact_comp.dualbound,...
%             info_mix1.dualbound, info_mix2.dualbound, info_mix3.dualbound,...
%             info_Linx.integrality_gap, info_DDFact.integrality_gap, info_DDFact_comp.integrality_gap,...
%             info_mix1.integrality_gap, info_mix2.integrality_gap, info_mix3.integrality_gap,...
%             info_Linx.continuous_dualgap, info_DDFact.continuous_dualgap, info_DDFact_comp.continuous_dualgap,...
%             info_mix1.continuous_dualgap, info_mix2.continuous_dualgap, info_mix3.continuous_dualgap,...
%             info_Linx.fixnum_to0, info_DDFact.fixnum_to0, info_DDFact_comp.fixnum_to0,...
%             info_mix1.fixnum_to0, info_mix2.fixnum_to0, info_mix3.fixnum_to0,...
%             info_Linx.fixnum_to1, info_DDFact.fixnum_to1, info_DDFact_comp.fixnum_to1,...
%             info_mix1.fixnum_to1, info_mix2.fixnum_to1, info_mix3.fixnum_to1,...
%             info_Linx.fixnum, info_DDFact.fixnum, info_DDFact_comp.fixnum,...
%             info_mix1.fixnum, info_mix2.fixnum, info_mix3.fixnum,...
%             info_mix1.alpha, info_mix2.alpha, info_mix3.alpha,...
%             info_Linx.time, info_DDFact.time, info_DDFact_comp.time,...
%             info_mix1.time, info_mix2.time, info_mix3.time,...
%             info_Linx.cputime, info_DDFact.cputime, info_DDFact_comp.cputime,...
%             info_mix1.cputime, info_mix2.cputime, info_mix3.cputime];
%     end
%     if repeattime==1
%         exceloutputall=exceloutput;
%     else
%         exceloutputall=exceloutputall+exceloutput;
%     end
% end
% exceloutputall=exceloutputall./totaltimes;
% title=["n", "s","lb",...
%     "info_Linx.exitflag", "info_DDFact.exitflag", "info_DDFact_comp.exitflag",...
%         "info_Linx.dualbound", "info_DDFact.dualbound", "info_DDFact_comp.dualbound",...
%         "info_mix1.dualbound", "info_mix2.dualbound", "info_mix3.dualbound",...
%         "info_Linx.integrality_gap", "info_DDFact.integrality_gap", "info_DDFact_comp.integrality_gap",...
%         "info_mix1.integrality_gap", "info_mix2.integrality_gap", "info_mix3.integrality_gap",...
%         "info_Linx.continuous_dualgap", "info_DDFact.continuous_dualgap", "info_DDFact_comp.continuous_dualgap",...
%         "info_mix1.continuous_dualgap", "info_mix2.continuous_dualgap", "info_mix3.continuous_dualgap",...
%         "info_Linx.fixnum_to0", "info_DDFact.fixnum_to0", "info_DDFact_comp.fixnum_to0",...
%         "info_mix1.fixnum_to0", "info_mix2.fixnum_to0", "info_mix3.fixnum_to0",...
%         "info_Linx.fixnum_to1", "info_DDFact.fixnum_to1", "info_DDFact_comp.fixnum_to1",...
%         "info_mix1.fixnum_to1", "info_mix2.fixnum_to1", "info_mix3.fixnum_to1",...
%         "info_Linx.fixnum", "info_DDFact.fixnum", "info_DDFact_comp.fixnum",...
%         "info_mix1.fixnum", "info_mix2.fixnum", "info_mix3.fixnum",...
%         "info_mix1.alpha", "info_mix2.alpha", "info_mix3.alpha",...
%         "info_Linx.time", "info_DDFact.time", "info_DDFact_comp.time",...
%         "info_mix1.time", "info_mix2.time", "info_mix3.time",...
%         "info_Linx.cputime", "info_DDFact.cputime", "info_DDFact_comp.cputime",...
%         "info_mix1.cputime", "info_mix2.cputime", "info_mix3.cputime"];
% xlswrite(fullFileNameexcel ,title,1,'A1');
% xlswrite(fullFileNameexcel ,exceloutputall,1,'A2');
% 
% 
