% %=============== Data63 =======================
% T63=readtable("./old_results/data63_mix_nonuniform_lincon_single.xlsx");
% T63H=readtable("./old_results/data63_single_lincon.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,47:48},T63{:,46}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% 
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("WallClock Time (seconds)");
% xlim([5,58]);
% % ylim([0,2.6])
% 
% legend("DDFact", "DDFactcomp","Linx","location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data63_nonuniform_lincon_time','-dpdf','-r0');
% 
% % %=============== Data63 =======================
% T63=readtable("./old_results/data63_mix_uniform_lincon_single.xlsx");
% T63H=readtable("./old_results/data63_single_lincon.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,47:48},T63{:,46}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% 
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("WallClock Time (seconds)");
% xlim([5,47]);
% % ylim([0,2.6])
% 
% legend("DDFact", "DDFactcomp","Linx","location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data63_uniform_lincon_time','-dpdf','-r0');
% %=============== Data63 END =======================

% %=============== Data63 =======================
% T63=readtable("./old_results/data63.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,47:48},T63{:,46}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% 
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("WallClock Time (seconds)");
% xlim([2,62]);
% % ylim([0,2.6])
% 
% legend("DDFact", "DDFactcomp","Linx","location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data63_time','-dpdf','-r0');

% %=============== Data63 END =======================


% %=============== Data90 =======================
% T63=readtable("./old_results/data90.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,47:48},T63{:,46}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% 
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("WallClock Time (seconds)");
% xlim([2,89]);
% % ylim([0,2.6])
% 
% legend("DDFact", "DDFactcomp","Linx","location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data90_time','-dpdf','-r0');
% %=============== Data90 END =======================
% % 
% %=============== Data124 =======================
% T63=readtable("./old_results/data124.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,47:48},T63{:,46}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% 
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("WallClock Time (seconds)");
% xlim([2,123]);
% % ylim([0,2.6])
% 
% legend("DDFact", "DDFactcomp","Linx","location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data124_time','-dpdf','-r0');
% %================ Data 124 end ========================
% % =================== Data 2000 ======================
T63=readtable("./old_results/data_mix2000_single.xlsx");
h63=figure;
h63.Position(3:4)=[1000,400];
nexttile;

p1=plot(T63{1:10,"s"},T63{1:10,[47,46]},"LineWidth",1.5);
p1(1).Color='#0072BD';
p1(2).Color='#EDB120';
xlabel("s");
ylabel("WallClock Time (seconds)");
xlim([20,200]);

nexttile;

p2=plot(T63{11:15,"s"},T63{11:15,[47,46]},"LineWidth",1.5);
p2(1).Color='#0072BD';
p2(2).Color='#EDB120';
xlabel("s");
ylabel("WallClock Time (seconds)");
xlim([860,940]);

legend("DDFact","Linx","location","northeastoutside");

set(h63,'Units','Inches');
pos = get(h63,'Position');
set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h63,'./graph/data2000_time','-dpdf','-r0');