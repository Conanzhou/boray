figure;
%%
subplots(2,1,1);hold on;
plot(yyray(1,1,1),yyray(1,3,1),'rx','linewidth',2);
contour(rr,zz,rhorz,'LevelList',0:0.05:0.7);hold on;
contour(rr,zz,rhorz,'LevelList',0.7:0.01:1.1);
for j=1:nray
    plot(yyray(:,1,j),yyray(:,3,j));
end
%%
num_x = length(yyray(:,1,1));
[point_rz_x,point_rz_o] = deal(ones(nray,1));
for j=2:nray
    [L,ind2]=deal(ones(num_x,1));
    for i=1000:num_x
        [L(i),ind2(i)]=leng([yyray(i,1,1),yyray(i,3,1)],[yyray(:,1,j),yyray(:,3,j)]);
    end
    [~,point_rz_x(j)] = min(L(:,1));
    point_rz_o(j) = ind2(point_rz_x(j));
    quiver(yyray(point_rz_x(j),1,1),yyray(point_rz_x(j),3,1),yyray(point_rz_x(j),4,1)*2e-5,yyray(point_rz_x(j),6,1)*2e-5,'b','LineWidth',2);
    quiver(yyray(point_rz_o(j),1,j),yyray(point_rz_o(j),3,j),yyray(point_rz_o(j),4,j)*2e-5,yyray(point_rz_o(j),6,j)*2e-5,'r','LineWidth',2);
end
ylim([0.0 0.03]);xlim([1.9 2.0])
TitleLabels('R-Z','','Z/m');set(gca,'xticklabel',[])
% run ./kB.m
%%
subplots(2,1,2);hold on;
tt=0:0.01*pi:2*pi;
plot(min(rg)*cos(tt),min(rg)*sin(tt),max(rg)*cos(tt),max(rg)*sin(tt));
for j=1:nray
    plot(yyray(:,1,j).*cos(yyray(:,2,j)),yyray(:,1,j).*sin(yyray(:,2,j)));
end
num_x = length(yyray(:,1,1));
[point_rphi_x,point_rphi_o] = deal(ones(nray,1));
for j=2:nray
    [L,ind2]=deal(ones(num_x,1));
    for i=1000:num_x
        [L(i),ind2(i)]=leng([yyray(i,1,1).*cos(yyray(i,2,1)),yyray(i,1,1).*sin(yyray(i,2,1))],[yyray(:,1,j).*cos(yyray(:,2,j)),yyray(:,1,j).*sin(yyray(:,2,j))]);
    end
    [~,point_rphi_x(j)] = min(L(:,1));
    point_rphi_o(j) = ind2(point_rphi_x(j));
    [x_temp1,y_temp1] = pol2cart(yyray(point_rphi_x(j),2,1) / 180 * pi,yyray(point_rphi_x(j),4,1));
    [x_temp2,y_temp2] = pol2cart((-1*yyray(point_rphi_x(j),2,1)+90) / 180 * pi,yyray(point_rphi_x(j),5,1)/yyray(point_rphi_x(j),1,1));
    quiver(yyray(point_rphi_x(j),1,1).*cos(yyray(point_rphi_x(j),2,1)),yyray(point_rphi_x(j),1,1).*sin(yyray(point_rphi_x(j),2,1)),(x_temp1+x_temp2)*2e-5,(y_temp1+y_temp2)*2e-5,'b','LineWidth',2);
    [x_temp1,y_temp1] = pol2cart(yyray(point_rphi_o(j),2,j) / 180 * pi,yyray(point_rphi_o(j),4,j));
    [x_temp2,y_temp2] = pol2cart((-1*yyray(point_rphi_o(j),2,j)+90) / 180 * pi,yyray(point_rphi_o(j),5,j)/yyray(point_rphi_o(j),1,j));
    quiver(yyray(point_rphi_o(j),1,j).*cos(yyray(point_rphi_o(j),2,j)),yyray(point_rphi_o(j),1,j).*sin(yyray(point_rphi_o(j),2,j)),(x_temp1+x_temp2)*2e-5,(y_temp1+y_temp2)*2e-5,'r','LineWidth',2);
end
ylim([0.0 0.07]);xlim([1.9 2.0])
TitleLabels('R-\Phi','R/m','\Phi/m');
% run ./kB_rphi.m

allAxes = findall(gcf,'type','axes');
linkaxes(allAxes,'x');
sgtitle(['f=',num2str(f/1e9,4),'GHz ']);
% set(gcf,'color','none');
% set(gcf,'InvertHardCopy','off');

%%
function [len,ind]=leng(p,l)
    L=sum( (l-p).^2 ,2);
    [len,ind]=min(L);
end
