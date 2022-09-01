% T63=readtable("./old_results/data63_mix_lincon.xlsx");
% T63H=readtable("./old_results/data63_single_lincon.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,38:39},T63{:,37},T63H{:,42}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% p(4).Color='#4DBEEE';
% p(4).LineStyle='--';
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("Number of Fixed Variables");
% xlim([2,58]);
% % ylim([0,3])
% 
% legend("DDFact", "DDFactcomp","Linx", "Mix DDFactcomp\_Linx", "location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data63_lincon_fixnum','-dpdf','-r0');

%===================================
% T63=readtable("./old_results/data63_mix_nonuniform_lincon_single.xlsx");
% T63H=readtable("./old_results/data63_single_lincon.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,38:39},T63{:,37}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("Number of Fixed Variables");
% xlim([5,58]);
% % ylim([0,3])
% 
% legend("DDFact", "DDFactcomp","Linx", "location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data63_nonuniform_lincon_fixnum','-dpdf','-r0');
% 
% %===================================
% T63=readtable("./old_results/data63_mix_uniform_lincon_single.xlsx");
% T63H=readtable("./old_results/data63_single_lincon.xlsx");
% 
% h63=figure;
% h63.Position(3:4)=[1000,400];
% p=plot(T63{:,"s"},[T63{:,38:39},T63{:,37}],"LineWidth",1.5);
% p(1).LineStyle="-";
% p(3).LineStyle="-";
% p(1).Color='#0072BD';
% p(2).Color='#D95319';
% p(3).Color='#EDB120';
% % title("Data63 Integrality Gap");
% xlabel("s");
% ylabel("Number of Fixed Variables");
% xlim([5,47]);
% % ylim([0,3])
% 
% legend("DDFact", "DDFactcomp","Linx", "location","northeastoutside");
% 
% set(h63,'Units','Inches');
% pos = get(h63,'Position');
% set(h63,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h63,'./graph/data63_uniform_lincon_fixnum','-dpdf','-r0');


% ============== Data 2000 ========
T2000=readtable("./old_results/data_mix2000_single.xlsx");

h2000=figure;
h2000.Position(3:4)=[1000,400];
nexttile;
p1=plot(T2000{1:10,"s"},T2000{1:10,[38]},"LineWidth",1.5);
p1(1).Color='#0072BD';
xlabel("s");
ylabel("Number of Fixed Variables");
xticks([20:20:200]);
xlim([20,200])
nexttile
p2=plot(T2000{11:15,"s"},T2000{11:15,[38]},"LineWidth",1.5);
p2(1).Color='#0072BD';
xlabel("s");
ylabel("Number of Fixed Variables");
ylim([0,10])

legend("DDFact","Linx","location","northeastoutside");
set(h2000,'Units','Inches');
pos = get(h2000,'Position');
set(h2000,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(h2000,'./graph/data2000_fixnum','-dpdf','-r0');
