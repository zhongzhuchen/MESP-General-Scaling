% load('data63.mat');
% % iterative_fixing_O;
% iterative_fixing_G;
% load('data90.mat');
% % iterative_fixing_O;
% iterative_fixing_G;
% load('data124.mat');
% % iterative_fixing_O;
% iterative_fixing_G;

%==========================================================================
load('data63.mat');
Linx_iterative_fixing_O;
load('data90.mat');
Linx_iterative_fixing_O;
load('data124.mat');
Linx_iterative_fixing_O;

%==========================================================================
% % solved_instance = struct;
% for n=[63, 90, 124]
%     filename = strcat('fix', num2str(n),'_LFcF_oscaling.mat');
%     load(filename);
%     % load('fix63_LFcF_oscaling.mat');
%     solvelist = [];
%     fixalllist1 = [];
%     fixnum1 = [];
%     for s=2:(n-1)
%         if fix(s).nearzero_intgap == 1
%             solvelist(end+1)=s;
%         end
%         if length(fix(s).fixto0list)+length(fix(s).fixto1list) >= n
%             fixalllist1(end+1)=s;
%         end
%         fixnum1(end+1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
%     end
%     fixalllist1 = union(solvelist, fixalllist1)
% %     solved_instance.oscaling_63_solvelist = fixalllist1;
% %     fixnums.oscaling_124 = fixnum;
% 
%     filename = strcat('fix', num2str(n),'_LFcF_gscaling.mat');
%     load(filename);
%     % load('fix63_LFcF_gscaling.mat');
%     solvelist = [];
%     fixalllist2 = [];
%     fixnum = [];
%     for s=2:(n-1)
%         if fix(s).nearzero_intgap == 1
%             solvelist(end+1)=s;
%         end
%         if length(fix(s).fixto0list)+length(fix(s).fixto1list) >= n
%             fixalllist2(end+1)=s;
%         end
%         fixnum(end+1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
%     end
%     fixalllist2 = union(solvelist, fixalllist2)
% %     solved_instance.gscaling_63_solvelist = fixalllist2;
% %     fixnums.gscaling_124 = fixnum;
% %     [c,~] = ismember(fixalllist1,fixalllist2)
% end
