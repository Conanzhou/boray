close all; 
clear; clc;

%% input wave parameters

% wave parameters
% k_phi0 = 180:2:200;

% f0 = ones(size(k_phi0))*(56); % GHz
f0 = 56:2:66;
num_f = size(f0,2); 
dt0_guess = ones(size(f0))*0.001; 
nt0_guess = ones(size(f0))*7000;

r = ones(size(f0))*2.2;
phih = ones(size(f0))*0;
z = ones(size(f0))*(-0.13);
% z = [-0.13 -0.14 -0.15 -0.16 -0.17];

% k_guess0 = [ -2000,-300];
k_theta0 = ones(size(f0))*(90);%度
k_phi0 = ones(size(f0))*180;

% equilibrium parameters
numeq = 1; % =0, analytical equilibrium; =1, numerical equilibrium
icase = 2;
eqfile = '../hl2a/input/genray_eqdata.mat'; % put it in the './input/' directory
% eqfile = '../eqdata/genray/EAST/genray_eqdata.mat'; % put it in the './input/' directory
savepath = '../hl2a/'; % '..' means from 'modules' directory

densprof(:,3)=0;
temprof(:,3)=0;


% [r,phih,z,kr_guess,nphi,kz,f]
% change kr_guess for possible multi-modes, i.e., O/X
yray0 = zeros(num_f,7);

for i = 1:num_f
    yray0(i,:) = [r(i),phih(i),z(i),-50,k_theta0(i),k_phi0(i),f0(i)];
%     yray0(i,:) = [2.2,0,-0.15,-600,nphi0(i),kz0(i),f0(i)];
%     yray0(i,:) = [2.2,0,-0.03,-1.845e+03,0,-2.97e+02,f0(i)];
end
%% calculate parameters

% dt0_guess = [0.001 0.001]; nt0_guess = [700 2000];
nray = size(yray0, 1);

%% calculate

disp('BO-Ray begins ray tracing.');
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
    k = fsolve(fDk, k_guess, options)
    kz = -1 * k * cos(k_theta/180*pi)
    kr = -1 *k * sin(k_theta/180*pi) * cos(k_phi/180*pi)
    nphi = -1 * k * sin(k_theta/180*pi) * sin(k_phi/180*pi) * r
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
    runtime1 = cputime;
    disp(['Begin to calculate the ray tracing use cold', ...
            ' plasma dispersion relation.']);
    run ./modules/boray_rk4;
%     yyray(:, :, jray) = yy;
    runtime1 = cputime - runtime1;
    disp(['Finish the ray tracing calculation. Runtime=', ...
            num2str(runtime1), 's']);
    
    %% plot the ray
    if(jray==1)
        fig_rphi=figure;hold on;
        tt=0:0.01*pi:2*pi;
        plot(min(rg)*cos(tt),min(rg)*sin(tt),max(rg)*cos(tt),max(rg)*sin(tt));
        axis equal;
        xlabel('X');ylabel('Y');

        title(['r=',num2str(yy(1,1),3),', \phi=',num2str(yy(1,2),3),', z=',num2str(yy(1,3),3)]);
        
        fig_rz=figure('unit','normalized','Position',[0.01 0.05 0.3 0.6],...
        'DefaultAxesFontSize',14);hold on;
        plot(yy(1,1).*cos(yy(1,2)),yy(1,1).*sin(yy(1,2)),'rx','linewidth',2);
        contour(rr,zz,fpsi,100);
        xlabel('R'); ylabel('Z');
        title(['f=',num2str(f/1e9,4),'GHz']);

    end
    figure(fig_rphi);
    plot(yy(:,1).*cos(yy(:,2)),yy(:,1).*sin(yy(:,2)));
    figure(fig_rz);
    plot(yy(:,1),yy(:,3),...
            yy(1,1),yy(1,3),'rx','linewidth',2);
end