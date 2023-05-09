% fixplot
% for n=[63, 90, 124]
%     y = zeros(1, n-2);
%     load(strcat('fix',num2str(63), '_LFcF.mat'))
%     for s = 2:(n-1)
%         y(s-1) = fix
%     end
% end

n = 63;
load('fix63_LFcF_gscaling.mat');
y1=zeros(n-2, 1);
y11=zeros(n-2,3);
for s = 2:(n-1)
    y1(s-1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
    y11(s-1,1) = min(length(fix(s).Linx_fixto0list_init)+length(fix(s).Linx_fixto1list_init), n);
    y11(s-1,2) = min(length(fix(s).Fact_fixto0list_init)+length(fix(s).Fact_fixto1list_init), n);
    y11(s-1,3) = min(length(fix(s).cFact_fixto0list_init)+length(fix(s).cFact_fixto1list_init), n);
end

load('fix63_LFcF_oscaling.mat');
y2=zeros(n-2, 1);
y22=zeros(n-2,3);
for s = 2:(n-1)
    y2(s-1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
    y22(s-1,1) = min(length(fix(s).Linx_fixto0list_init)+length(fix(s).Linx_fixto1list_init), n);
    y22(s-1,2) = min(length(fix(s).Fact_fixto0list_init)+length(fix(s).Fact_fixto1list_init), n);
    y22(s-1,3) = min(length(fix(s).cFact_fixto0list_init)+length(fix(s).cFact_fixto1list_init), n);
end

y=[y1, y2];
h=figure;
h.Position(3:4)=[1000,400];
p = plot(20:39, y(20:39,:), "LineWidth",1.5);
p(1).Color ="#000000";
p(2).Color ="#000000";
p(2).LineStyle = "--";
lgnd = legend("G-scale", "O-scale", "location","northeast");
set(lgnd,'color','none');
title(strcat("Number of fixed variables: G-scaling vs O-scaling (n=", num2str(n),")"));
xlabel("s");
ylabel("Number of fixed variables");
xticks([20, 30, 40]);
xlim([20,43]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('./graph/data', num2str(n), 'fix'),'-dpdf','-r0');

% y=[y11, y22];
% h=figure;
% h.Position(3:4)=[1000,400];
% p = plot(2:(n-1), y, "LineWidth",1.5);
% p(1).Color ="#EDB120";
% p(2).Color ="#0072BD";
% p(3).Color ="#D95319";
% p(4).Color ="#EDB120";
% p(5).Color ="#0072BD";
% p(6).Color ="#D95319";
% for i = 4:6
%     p(i).LineStyle = "--";
% end
% lgnd = legend("G-scale-Linx", "G-scale-Fact", "G-scale-cFact",...
%     "O-scale-Linx", "O-scale-Fact", "O-scale-cFact", "location","northeast");
% set(lgnd,'color','none');
% title(strcat("Number of fixed variables: G-scaling vs O-scaling (n=", num2str(n),")"));
% xlabel("s");
% ylabel("Number of fixed variables");
% xticks([0, 10, 20, 30, 40, 50, 60]);
% xlim([0,73]);
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,strcat('./graph/data', num2str(n), 'separate-fix'),'-dpdf','-r0');

%--------------------------------------------------------------------------
n = 90;
load('fix90_LFcF_gscaling.mat');
y1=zeros(n-2, 1);
y11=zeros(n-2,3);
for s = 2:(n-1)
    y1(s-1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
    y11(s-1,1) = min(length(fix(s).Linx_fixto0list_init)+length(fix(s).Linx_fixto1list_init), n);
    y11(s-1,2) = min(length(fix(s).Fact_fixto0list_init)+length(fix(s).Fact_fixto1list_init), n);
    y11(s-1,3) = min(length(fix(s).cFact_fixto0list_init)+length(fix(s).cFact_fixto1list_init), n);
end

load('fix90_LFcF_oscaling.mat');
y2=zeros(n-2, 1);
y22=zeros(n-2,3);
for s = 2:(n-1)
    y2(s-1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
    y22(s-1,1) = min(length(fix(s).Linx_fixto0list_init)+length(fix(s).Linx_fixto1list_init), n);
    y22(s-1,2) = min(length(fix(s).Fact_fixto0list_init)+length(fix(s).Fact_fixto1list_init), n);
    y22(s-1,3) = min(length(fix(s).cFact_fixto0list_init)+length(fix(s).cFact_fixto1list_init), n);
end

y=[y1, y2];
h=figure;
h.Position(3:4)=[1000,400];
p = plot(22:68, y(22:68,:), "LineWidth",1.5);
p(1).Color ="#000000";
p(2).Color ="#000000";
p(2).LineStyle = "--";
lgnd = legend("G-scale", "O-scale", "location","northeast");
set(lgnd,'color','none');
title(strcat("Number of fixed variables: G-scaling vs O-scaling (n=", num2str(n),")"));
xlabel("s");
ylabel("Number of fixed variables");
xticks([20, 30, 40, 50, 60, 70]);
xlim([20,78]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('./graph/data', num2str(n), 'fix'),'-dpdf','-r0');

% y=[y11, y22];
% h=figure;
% h.Position(3:4)=[1000,400];
% p = plot(2:(n-1), y, "LineWidth",1.5);
% p(1).Color ="#EDB120";
% p(2).Color ="#0072BD";
% p(3).Color ="#D95319";
% p(4).Color ="#EDB120";
% p(5).Color ="#0072BD";
% p(6).Color ="#D95319";
% for i = 4:6
%     p(i).LineStyle = "--";
% end
% lgnd = legend("G-scale-Linx", "G-scale-Fact", "G-scale-cFact",...
%     "O-scale-Linx", "O-scale-Fact", "O-scale-cFact", "location","northeast");
% set(lgnd,'color','none');
% title(strcat("Number of fixed variables: G-scaling vs O-scaling (n=", num2str(n),")"));
% xlabel("s");
% ylabel("Number of fixed variables");
% xticks([0, 10, 20, 30, 40, 50, 60, 70,80, 90]);
% xlim([0,100]);
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('./graph/data', num2str(n), 'separate-fix'),'-dpdf','-r0');

%--------------------------------------------------------------------------
n = 124;
load('fix124_LFcF_gscaling.mat');
y1=zeros(n-2, 1);
y11=zeros(n-2,3);
for s = 2:(n-1)
    y1(s-1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
    y11(s-1,1) = min(length(fix(s).Linx_fixto0list_init)+length(fix(s).Linx_fixto1list_init), n);
    y11(s-1,2) = min(length(fix(s).Fact_fixto0list_init)+length(fix(s).Fact_fixto1list_init), n);
    y11(s-1,3) = min(length(fix(s).cFact_fixto0list_init)+length(fix(s).cFact_fixto1list_init), n);
end

load('fix124_LFcF_oscaling.mat');
y2=zeros(n-2, 1);
y22=zeros(n-2,3);
for s = 2:(n-1)
    y2(s-1) = min(length(fix(s).fixto0list)+length(fix(s).fixto1list), n);
    y22(s-1,1) = min(length(fix(s).Linx_fixto0list_init)+length(fix(s).Linx_fixto1list_init), n);
    y22(s-1,2) = min(length(fix(s).Fact_fixto0list_init)+length(fix(s).Fact_fixto1list_init), n);
    y22(s-1,3) = min(length(fix(s).cFact_fixto0list_init)+length(fix(s).cFact_fixto1list_init), n);
end

y=[y1, y2];
h=figure;
h.Position(3:4)=[1000,400];
p = plot(27:113, y(27:113,:), "LineWidth",1.5);
p(1).Color ="#000000";
p(2).Color ="#000000";
p(2).LineStyle = "--";
lgnd = legend("G-scale", "O-scale", "location","northeast");
set(lgnd,'color','none');
title(strcat("Number of fixed variables: G-scaling vs O-scaling (n=", num2str(n),")"));
xlabel("s");
ylabel("Number of fixed variables");
xticks([20, 30, 40, 50, 60, 70,80, 90, 100, 110, 120]);
xlim([27,135]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('./graph/data', num2str(n), 'fix'),'-dpdf','-r0');

% y=[y11, y22];
% h=figure;
% h.Position(3:4)=[1000,400];
% p = plot(2:(n-1), y, "LineWidth",1.5);
% p(1).Color ="#EDB120";
% p(2).Color ="#0072BD";
% p(3).Color ="#D95319";
% p(4).Color ="#EDB120";
% p(5).Color ="#0072BD";
% p(6).Color ="#D95319";
% for i = 4:6
%     p(i).LineStyle = "--";
% end
% lgnd = legend("G-scale-Linx", "G-scale-Fact", "G-scale-cFact",...
%     "O-scale-Linx", "O-scale-Fact", "O-scale-cFact", "location","northeast");
% set(lgnd,'color','none');
% title(strcat("Number of fixed variables: G-scaling vs O-scaling (n=", num2str(n),")"));
% xlabel("s");
% ylabel("Number of fixed variables");
% xticks([0, 10, 20, 30, 40, 50, 60, 70,80, 90, 100, 110, 120]);
% xlim([0,145]);
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,strcat('./graph/data', num2str(n), 'separate-fix'),'-dpdf','-r0');


