% load('yyray.mat', 'yyray');
% figure();hold on;
%%
xmode1 = yyray(1:18000,:,1);
x1x = xmode1(:,1);
x1y = xmode1(:,3);
% plot(x1x,x1y,'.');
%%
omode1 = yyray(1:6000,:,2);
o1x = omode1(:,1);
o1y = omode1(:,3);
% plot(o1x,o1y,'.');
%
omode2 = yyray(1:6000,:,3);
o2x = omode2(:,1);
o2y = omode2(:,3);
% plot(o2x,o2y,'.');
%
omode3 = yyray(1:6000,:,4);
o3x = omode3(:,1);
o3y = omode3(:,3);
% plot(o3x,o3y,'.');
%%
%插值法

%%
%拟合法,y-x
facter0 = polyfit(x1y,x1x,6);
facter1 = polyfit(o1y,o1x,3);
facter2 = polyfit(o2y,o2x,3);
facter3 = polyfit(o3y,o3x,3);

syms x y
sol = solve(y>0.005,x==facter1(1)*y^3+facter1(2)*y^2+facter1(3)*y^1+facter1(4),x==facter0(1)*y^6+facter0(2)*y^5+facter0(3)*y^4+facter0(4)*y^3+facter0(5)*y^2+facter0(6)*y^1+facter0(7));
cross1_x = double(sol.x);
cross1_y = double(sol.y);
%
sol = solve(y>0.005,x==facter2(1)*y^3+facter2(2)*y^2+facter2(3)*y^1+facter2(4),x==facter0(1)*y^6+facter0(2)*y^5+facter0(3)*y^4+facter0(4)*y^3+facter0(5)*y^2+facter0(6)*y^1+facter0(7));
cross2_x = double(sol.x);
cross2_y = double(sol.y);
%
sol = solve(y>0.005,x==facter3(1)*y^3+facter3(2)*y^2+facter3(3)*y^1+facter3(4),x==facter0(1)*y^6+facter0(2)*y^5+facter0(3)*y^4+facter0(4)*y^3+facter0(5)*y^2+facter0(6)*y^1+facter0(7));
cross3_x = double(sol.x);
cross3_y = double(sol.y);

%% 
o1_num = find(abs(o1y-cross1_y)<=1e-6);
x1_num = find(abs(x1y-cross1_y)<=0.3e-6);
% find(abs(o1y-cross_y)<=0.00001)
x1_cross_ray = yyray(x1_num(1),:,1);
o1_cross_ray = yyray(o1_num(1),:,2);
quiver(o1_cross_ray(1),o1_cross_ray(3),o1_cross_ray(4)*2e-5,o1_cross_ray(6)*2e-5,'b','LineWidth',2);
quiver(x1_cross_ray(1),x1_cross_ray(3),x1_cross_ray(4)*2e-5,x1_cross_ray(6)*2e-5,'r','LineWidth',2);
%
o2_num = find(abs(o2y-cross2_y)<=1.2e-6);
x2_num = find(abs(x1y-cross2_y)<=0.3e-6);
% find(abs(o2y-cross_y)<=0.00001)
x2_cross_ray = yyray(x2_num(1),:,1);
o2_cross_ray = yyray(o2_num(1),:,3);
quiver(o2_cross_ray(1),o2_cross_ray(3),o2_cross_ray(4)*2e-5,o2_cross_ray(6)*2e-5,'b','LineWidth',2);
quiver(x2_cross_ray(1),x2_cross_ray(3),x2_cross_ray(4)*2e-5,x2_cross_ray(6)*2e-5,'r','LineWidth',2);
%
o3_num = find(abs(o3y-cross3_y)<=3e-6);
x3_num = find(abs(x1y-cross3_y)<=1.5e-6);
% find(abs(o3y-cross_y)<=0.00001)
x3_cross_ray = yyray(x3_num(1),:,1);
o3_cross_ray = yyray(o3_num(1),:,4);
quiver(o3_cross_ray(1),o3_cross_ray(3),o3_cross_ray(4)*2e-5,o3_cross_ray(6)*2e-5,'b','LineWidth',2);
quiver(x3_cross_ray(1),x3_cross_ray(3),x3_cross_ray(4)*2e-5,x3_cross_ray(6)*2e-5,'r','LineWidth',2);