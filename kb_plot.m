figure;
% subplots(2,1,1);
subplots(1,3,1);
set(gca,'DataAspectRatio',[1 1 1])
hold on;
plot(yyray(1,1,1),yyray(1,3,1),'rx','linewidth',2);
contour(rr,zz,rhorz,'LevelList',0:0.05:0.7);hold on;
contour(rr,zz,rhorz,'LevelList',0.7:0.01:1.1);
plot(yyray(:,1,1),yyray(:,3,1),'b','linewidth',2);

%%
num_x = length(yyray(:,1,1));
[point_rz_x,point_rz_o] = deal(ones(nray,1));
rz_output_rho = ones(nray-1,2);
rz_output_k_r = ones(nray-1,2); 
rz_output_k_phi = ones(nray-1,2); 
rz_output_k_z = ones(nray-1,2); 

for j=2:nray
    [L,ind2]=deal(ones(num_x,1));
    for i=1000:num_x
        [L(i),ind2(i)]=leng([yyray(i,1,1),yyray(i,3,1)],[yyray(:,1,j),yyray(:,3,j)]);
    end
    [~,point_rz_x(j)] = min(L(:,1));
    point_rz_o(j) = ind2(point_rz_x(j));
    
    plot(yyray(1:point_rz_o(j),1,j),yyray(1:point_rz_o(j),3,j),'r','linewidth',2);
    
    rz_output_rho(j-1,1) = interp2(rr',zz',rhorz',yyray(point_rz_x(j),1,1),yyray(point_rz_x(j),3,1));
    rz_output_rho(j-1,2) = interp2(rr',zz',rhorz',yyray(point_rz_o(j),1,j),yyray(point_rz_o(j),3,j));
    rz_output_k_r(j-1,1) = yyray(point_rz_x(j),4,1);
    rz_output_k_r(j-1,2) = -1*yyray(point_rz_o(j),4,j);
    rz_output_k_phi(j-1,1) = yyray(point_rz_x(j),5,1) / yyray(point_rz_x(j),1,1);
    rz_output_k_phi(j-1,2) = -1*yyray(point_rz_o(j),5,j) / yyray(point_rz_o(j),1,j);
    rz_output_k_z(j-1,1) = yyray(point_rz_x(j),6,1);
    rz_output_k_z(j-1,2) = -1*yyray(point_rz_o(j),6,j);

    quiver(yyray(point_rz_x(j),1,1),yyray(point_rz_x(j),3,1),yyray(point_rz_x(j),4,1)*2e-5,yyray(point_rz_x(j),6,1)*2e-5,'b','LineWidth',3);
    quiver(yyray(point_rz_o(j),1,j),yyray(point_rz_o(j),3,j),-1*yyray(point_rz_o(j),4,j)*2e-5,-1*yyray(point_rz_o(j),6,j)*2e-5,'r','LineWidth',3);
end

% ylim([0.0 0.03]);xlim([1.9 2.0])
% ylim([0.0 0.04]);xlim([2.3 2.45])
TitleLabels('(a)CPS,k_B = k_O-k_X','R/m','Z/m');%set(gca,'xticklabel',[])
% run ./kB.m
%%
rz_kb_rho = (rz_output_rho(:,1)+rz_output_rho(:,2))/2;
rz_kb_r = rz_output_k_r(:,2) - rz_output_k_r(:,1);
rz_kb_phi = rz_output_k_phi(:,2) - rz_output_k_phi(:,1);
rz_kb_z = rz_output_k_z(:,2) - rz_output_k_z(:,1);
%%
subplots(2,1,2);hold on;
tt=0:0.01*pi:2*pi;
plot(min(rg)*cos(tt),min(rg)*sin(tt),max(rg)*cos(tt),max(rg)*sin(tt));
plot(yyray(1,1,1).*cos(yyray(1,2,1)),yyray(1,1,1).*sin(yyray(1,2,1)),'rx','linewidth',2);
for j=1:nray
    plot(yyray(:,1,j).*cos(yyray(:,2,j)),yyray(:,1,j).*sin(yyray(:,2,j)));
end
num_x = length(yyray(:,1,1));
[point_rphi_x,point_rphi_o] = deal(ones(nray,1));
rphi_output_rho = ones(nray-1,2);
rphi_output_k_r = ones(nray-1,2); 
rphi_output_k_phi = ones(nray-1,2); 
rphi_output_k_z = ones(nray-1,2); 

for j=2:nray
    [L,ind2]=deal(ones(num_x,1));
    for i=1000:num_x
        [L(i),ind2(i)]=leng([yyray(i,1,1).*cos(yyray(i,2,1)),yyray(i,1,1).*sin(yyray(i,2,1))],[yyray(:,1,j).*cos(yyray(:,2,j)),yyray(:,1,j).*sin(yyray(:,2,j))]);
    end
    [~,point_rphi_x(j)] = min(L(:,1));
    point_rphi_o(j) = ind2(point_rphi_x(j));

    rphi_output_rho(j-1,1) = interp2(rr',zz',rhorz',yyray(point_rphi_x(j),1,1),yyray(point_rphi_x(j),3,1));
    rphi_output_rho(j-1,2) = interp2(rr',zz',rhorz',yyray(point_rphi_o(j),1,j),yyray(point_rphi_o(j),3,j));
    rphi_output_k_r(j-1,1) = yyray(point_rphi_x(j),4,1);
    rphi_output_k_r(j-1,2) = -1*yyray(point_rphi_o(j),4,j);
    rphi_output_k_phi(j-1,1) = yyray(point_rphi_x(j),5,1) / yyray(point_rphi_x(j),1,1);
    rphi_output_k_phi(j-1,2) = -1*yyray(point_rphi_o(j),5,j) / yyray(point_rphi_o(j),1,j);
    rphi_output_k_z(j-1,1) = yyray(point_rphi_x(j),6,1);
    rphi_output_k_z(j-1,2) = -1*yyray(point_rphi_o(j),6,j);


    [x_temp1,y_temp1] = pol2cart(yyray(point_rphi_x(j),2,1),yyray(point_rphi_x(j),4,1));
    [x_temp2,y_temp2] = pol2cart(yyray(point_rphi_x(j),2,1)+pi/2,yyray(point_rphi_x(j),5,1)/yyray(point_rphi_x(j),1,1));
    quiver(yyray(point_rphi_x(j),1,1).*cos(yyray(point_rphi_x(j),2,1)),yyray(point_rphi_x(j),1,1).*sin(yyray(point_rphi_x(j),2,1)),(x_temp1+x_temp2)*2e-5,(y_temp1+y_temp2)*2e-5,'b','LineWidth',2);
    [x_temp1,y_temp1] = pol2cart(yyray(point_rphi_o(j),2,j),yyray(point_rphi_o(j),4,j));
    [x_temp2,y_temp2] = pol2cart(yyray(point_rphi_o(j),2,j)+pi/2,yyray(point_rphi_o(j),5,j)/yyray(point_rphi_o(j),1,j));
    quiver(yyray(point_rphi_o(j),1,j).*cos(yyray(point_rphi_o(j),2,j)),yyray(point_rphi_o(j),1,j).*sin(yyray(point_rphi_o(j),2,j)),-1*(x_temp1+x_temp2)*2e-5,-1*(y_temp1+y_temp2)*2e-5,'r','LineWidth',2);
end
% ylim([0.0 0.07]);xlim([1.9 2.0])
% ylim([0.0 0.04]);xlim([2.3 2.45])
TitleLabels('R-\Phi','X/m','Y/m');
% run ./kB_rphi.m
%%
rphi_kb_rho = (rphi_output_rho(:,1)+rphi_output_rho(:,2))/2;
rphi_kb_r = rphi_output_k_r(:,2) - rphi_output_k_r(:,1);
rphi_kb_phi = rphi_output_k_phi(:,2) - rphi_output_k_phi(:,1);
rphi_kb_z = rphi_output_k_z(:,2) - rphi_output_k_z(:,1);
%%
allAxes = findall(gcf,'type','axes');
linkaxes(allAxes,'x');
sgtitle(['f=',num2str(f/1e9,4),'GHz ']);
% set(gcf,'color','none');
% set(gcf,'InvertHardCopy','off');
%%
% figure;
open("./hl2m/output/kB_spectrum.fig");
hold on;
rz_kb = sqrt(rz_kb_r.^2+rz_kb_phi.^2+rz_kb_z.^2);
plot(rz_kb_rho,rz_kb/100,'-o','linewidth',2,'MarkerSize',5);
TitleLabels 'kB spectrum' '\rho/a' 'k_B/cm^{-1}'
%%
function [len,ind]=leng(p,l)
    L=sum( (l-p).^2 ,2);
    [len,ind]=min(L);
end
