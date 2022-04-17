close all; 
clear; clc;


% input wave parameters
f0 = 50:2:65;
% f0 = 20:2.5:35;
% f0 = 51;
% k_theta0 = 35:5:130;%度,范围选取，基于平衡位型，基本相切,down
% k_theta0 = 40:5:145;%度,范围选取，基于平衡位型，基本相切,up
% k_theta0 = 35:5:145;%度,范围选取，基于平衡位型，基本相切,中平面
% k_theta0 = 90;
k_theta0 = 60:5:120;

num_f = size(f0,2); 
num_theta = size(k_theta0,2); 

dt0_guess = ones(num_f*num_theta)*0.001;
nt0_guess = ones(num_f*num_theta)*2100;
% nt0_guess = ones(num_f*num_theta)*2100*2;

r = ones(num_f*num_theta)*2.63842;
phih = ones(num_f*num_theta)*0;
% z = ones(num_f*num_theta)*(0.158);
z = ones(num_f*num_theta)*(-0.158);
% z = [ones([1 num_f/2])*(-0.158) ones([1 num_f/2])*(0.158)];
% z = [0 0.001 0.0015 0.002];

% k_guess0 = [ -500,-2000,-2000,-2000];
% k_guess0 = [ -500 -2000*ones(1,num_f)];
k_guess0 = ones(num_f*num_theta) * (-10);
% k_guess0 = ones(num_f*num_theta) * (-2000);
% k_guess0 = [ones([1 num_f/2])*(-10) ones([1 num_f/2])*(-2000)];
% k_theta0 = ones(size(f0))*(35);%度
% k_theta0 = [88 87 86 85]-5;
% k_theta0 = [84 83.9 83.5 82]-4;
% k_theta0 = [84 83.5 83.5 83.5]+4;
% k_theta0 = [ones([1 num_f/2])*(100) ones([1 num_f/2])*(135)];
k_phi0 = ones(num_f*num_theta)*180;
% k_phi0 = ones(size(f0))*175;
% k_phi0 = [175 172 171  169];

% equilibrium parameters
numeq = 1; % =0, analytical equilibrium; =1, numerical equilibrium
icase = 2;
% eqfile = '../hl2m/input/L_mode/genray_eqdata.mat'; % put it in the './input/' directory
eqfile = '../hl2m/input/genray_eqdata.mat'; % put it in the './input/' directory
% eqfile = '../eqdata/genray/EAST/genray_eqdata.mat'; % put it in the './input/' directory
savepath = './hl2m/output/H_mode/X-mode/down/'; % '..' means from 'modules' directory

densprof(:,3)=0;
temprof(:,3)=0;


% [r,phih,z,kr_guess,nphi,kz,f]
% change kr_guess for possible multi-modes, i.e., O/X
yray0 = zeros(num_f*num_theta,7);

for i = 1:num_f
    for j = 1:num_theta
        temp_num = (i-1)*num_theta+j;
        yray0(temp_num,:) = [r(temp_num),phih(temp_num),z(temp_num),k_guess0(temp_num),k_theta0(j),k_phi0(temp_num),f0(i)];
%       yray0(i,:) = [2.2,0,-0.15,-600,nphi0(i),kz0(i),f0(i)];
%       yray0(i,:) = [2.2,0,-0.03,-1.845e+03,0,-2.97e+02,f0(i)];
    end
end
%% calculate parameters

% dt0_guess = [0.001 0.001]; nt0_guess = [700 2000];a
nray = size(yray0, 1);

%% calculate

disp('BO-Ray begins ray tracing.');
runtime1 = cputime;
for jray = 1:nray
    %%
    r = yray0(jray, 1); phi = yray0(jray, 2); z = yray0(jray, 3);
    k_guess = yray0(jray, 4); k_theta = yray0(jray, 5); k_phi = yray0(jray, 6);
    f = yray0(jray,7)*1e9;
    run ./modules/boray_initialize2;
    disp('Finish the initialize. Later to solve the initial kr.');

    %%
    % # multi-solutions may exists, choice the correct initial value to
    % determine which mode to tracing
    options = optimoptions('fsolve', 'Display', 'none', 'TolX', 1e-8);
    % kr=fsolve(fDkr,-500.0,options) % O-mode
    % kr=fsolve(fDkr,-0.1,options) % X-mode
    k = fsolve(fDk, k_guess, options);
    kz = -1 * k * cos(k_theta/180*pi);
    kr = -1 *k * sin(k_theta/180*pi) * cos(k_phi/180*pi);
    nphi = -1 * k * sin(k_theta/180*pi) * sin(k_phi/180*pi) * r;
    %%
    % # calculate the ray tracing use RK4
    dt0 = dt0_guess(jray);
    nt0 = nt0_guess(jray);
    tmp = 20; % rescaling the dt
    nt = nt0 * tmp; dt = dt0 / c / tmp; % time steps and total time
    if (jray == 1)
        yyray = zeros(nt, 18 + 2 * S, nray) + NaN; %  store for ray tracing
        yrayp = zeros(nt, 6 + (1 + S), nray) + NaN; % store for power absorb
    end
%     runtime1 = cputime;
%     disp(['Begin to calculate the ray tracing use cold', ...
%             ' plasma dispersion relation.']);
    disp(['ray=',num2str(jray)]);
    run ./modules/boray_rk4_SZH;
    yyray(:, :, jray) = yy;
%     runtime1 = cputime - runtime1;
%     disp(['Finish the ray tracing calculation. Runtime=', ...
%             num2str(runtime1), 's']);
    
%     %% plot the ray
%     if(jray==1)
%         fig_rphi=figure;hold on;
%         tt=0:0.01*pi:2*pi;
%         plot(min(rg)*cos(tt),min(rg)*sin(tt),max(rg)*cos(tt),max(rg)*sin(tt));
%         axis equal;box on;
%         xlabel('X');ylabel('Y');
% 
%         title(['r=',num2str(yy(1,1),3),', \phi=',num2str(yy(1,2),3),', z=',num2str(yy(1,3),3)]);
%         
%         plot(yy(1,1).*cos(yy(1,2)),yy(1,1).*sin(yy(1,2)),'rx','linewidth',2);
%         fig_rz=figure('unit','normalized','Position',[0.01 0.05 0.3 0.6],'DefaultAxesFontSize',14);hold on;
%         ylim([-1.2 1.2]);
% %         contour(rr,zz,fpsi,100);
%         plot(yy(1,1),yy(1,3),'rx','linewidth',2);hold on;
% %         contour(rr,zz,rhorz,'LevelList',0:0.05:0.7);hold on;
% %         contour(rr,zz,rhorz,'LevelList',0.7:0.01:1.1);
% %         contour(rr,zz,squeeze(fns0(1,:,:)),100);
%         contour(rr,zz,squeeze(fns0(1,:,:)),'LevelList',(0.6:0.01:2)*1e19);
%         xlabel('R'); ylabel('Z');box on;
%         title(['f=',num2str(f/1e9,4),'GHz ','theta=',num2str(k_theta0(1))]);
%         title(['f=52:1:65GHz ','theta=',num2str(k_theta0(1)-90)]);
% 
%     end
%     figure(fig_rphi);
%     plot(yy(:,1).*cos(yy(:,2)),yy(:,1).*sin(yy(:,2)));%ylim([0.0 0.07]);xlim([1.9 2.05]);
%     figure(fig_rz);
%     plot(yy(:,1),yy(:,3),'.');
end
%%
    runtime1 = cputime - runtime1;
    disp(['Finish the ray tracing calculation. Runtime=', ...
            num2str(runtime1), 's']);
% run ./kB.m
% ylim([0.0 0.07]);xlim([1.9 2.05])
% figure(fig_rphi);
% run ./kB_rphi.m
% ylim([0.0 0.07]);xlim([1.9 2.05]);
%%
% run ./k_spectrum.m
% run ./kper_plot_HTPD.m