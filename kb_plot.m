figure;
subplots(2,1,1);hold on;
plot(yyray(1,1,1),yyray(1,3,1),'rx','linewidth',2);
contour(rr,zz,rhorz,'LevelList',0:0.05:0.7);hold on;
contour(rr,zz,rhorz,'LevelList',0.7:0.01:1.1);
for j=1:nray
    plot(yyray(:,1,j),yyray(:,3,j));
end
ylim([0.0 0.03]);xlim([1.9 2.0])
TitleLabels('R-Z','','Z/m');set(gca,'xticklabel',[])
run ./kB.m

subplots(2,1,2);hold on;
tt=0:0.01*pi:2*pi;
plot(min(rg)*cos(tt),min(rg)*sin(tt),max(rg)*cos(tt),max(rg)*sin(tt));
for j=1:nray
    plot(yyray(:,1,j).*cos(yyray(:,2,j)),yyray(:,1,j).*sin(yyray(:,2,j)));
end
ylim([0.0 0.07]);xlim([1.9 2.0])
TitleLabels('R-\Phi','R/m','\Phi/m');
run ./kB_rphi.m

allAxes = findall(gcf,'type','axes');
linkaxes(allAxes,'x');
sgtitle(['f=',num2str(f/1e9,4),'GHz ']);
set(gcf,'color','none');
set(gcf,'InvertHardCopy','off');