% close all; 
clear; clc;

%% equilibrium parameters

numeq = 1; % =0, analytical equilibrium; =1, numerical equilibrium
eqfile = '../hl2a/input/genray_eqdata.mat'; % put it in the './input/' directory
densprof(:,3)=0;
temprof(:,3)=0;
%% wave parameters

savepath = '../hl2a/output'; % '..' means from 'modules' directory

densprof(:,3)=0;
temprof(:,3)=0;
%% wave parameters

f0 = 20:0.5:35; % GHz
% k_guess0 = [ -2000,-300];
num_f = size(f0,2); 
r00 = ones(size(f0))*2.2;
phih00 = ones(size(f0))*0;
z00 = ones(size(f0))*(0);
% z = [-0.13 -0.14 -0.15 -0.16 -0.17];

k_theta0 = ones(size(f0))*(90);%度
k_phi0 = ones(size(f0))*(180);

% [r,phih,z,kr_guess,nphi,kz,f]
% change kr_guess for possible multi-modes, i.e., O/X
yray0 = zeros(num_f,7);

for i = 1:num_f
    yray0(i,:) = [r00(i),phih00(i),z00(i),-10,k_theta0(i),k_phi0(i),f0(i)];
end
%% calculate parameters
nray = size(yray0, 1);
%% calculate
disp('BO-Ray begins ray tracing.');
for jray = 1:nray
    %%54
    r = yray0(jray, 1); phi = yray0(jray, 2); z = yray0(jray, 3);
    k_guess = yray0(jray, 4); k_theta = yray0(jray, 5); k_phi = yray0(jray, 6);
    f = yray0(jray,7)*1e9;
    run ./modules/boray_initialize2;
    disp('Finish the initialize. Later to solve the initial kr.');
    %%
    % # multi-solutions may exists, choice the correct initial value to
    % determine which mode to tracing
    options = optimoptions('fsolve', 'Display', 'none', 'TolX', 1e-9);
    % kr=fsolve(fDkr,-500.0,options) % O-mode
    % kr=fsolve(fDkr,-0.1,options) % X-mode
    k = fsolve(fDk, k_guess, options)
    ans_k(jray) = -k;
end
%%
figure;hold on;
plot(ans_k,f0,'ro');

%%
for i = 1:num_f
    yray0(i,:) = [r00(i),phih00(i),z00(i),-2000,k_theta0(i),k_phi0(i),f0(i)];
end
nray = size(yray0, 1);
for jray = 1:nray
    r = yray0(jray, 1); phi = yray0(jray, 2); z = yray0(jray, 3);
    k_guess = yray0(jray, 4); k_theta = yray0(jray, 5); k_phi = yray0(jray, 6);
    f = yray0(jray,7)*1e9;
    run ./modules/boray_initialize2;
    options = optimoptions('fsolve', 'Display', 'none', 'TolX', 1e-8);
    k = fsolve(fDk, k_guess, options)
    ans_k(jray) = -k;
end
plot(ans_k,f0,'b*');xlabel('k'); ylabel('f');