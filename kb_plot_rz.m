% figure;
% subplots(2,1,1);
subplots(1,3,1);
set(gca,'DataAspectRatio',[1 1 1])
% hold on;
plot(yyray(1,1,1),yyray(1,3,1),'rx','linewidth',3);
contour(rr,zz,rhorz,'LevelList',0:0.05:0.7);hold on;
contour(rr,zz,rhorz,'LevelList',0.7:0.01:1.1);
plot(yyray(:,1,1),yyray(:,3,1),'b','linewidth',3);

for j=2:nray
    
    plot(yyray(1:point_rz_o(j),1,j),yyray(1:point_rz_o(j),3,j),'r','linewidth',3);
    
    quiver(yyray(point_rz_x(j),1,1),yyray(point_rz_x(j),3,1),yyray(point_rz_x(j),4,1)*2e-5,yyray(point_rz_x(j),6,1)*2e-5,'g','LineWidth',1,'MaxHeadSize',5,'AutoScaleFactor',2);
    quiver(yyray(point_rz_o(j),1,j),yyray(point_rz_o(j),3,j),-1*yyray(point_rz_o(j),4,j)*2e-5,-1*yyray(point_rz_o(j),6,j)*2e-5,'y','LineWidth',1,'MaxHeadSize',5,'AutoScaleFactor',2);
end
