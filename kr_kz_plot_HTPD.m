%%
% load('D:\0CODE\boray\hl2m\input\genray_eqdata.mat')
% load('D:\0CODE\boray\hl2m\output\up\k_output_DBS_H_UP.mat')
% load('D:\Temp\0data\HTPD\H-mode\X-mode\yyray_up.mat')
f0 = 50:2:65;
% k_theta0 = 35:5:145;
k_theta0 = 60:5:90;
num_f = size(f0,2); 
num_theta = size(k_theta0,2); 
%%
k_output = zeros(num_f,num_theta,8);
for i=1:num_f
    for j = 1:num_theta
        temp_num = (i-1)*num_theta+j;
        [k_output(i,j,1),ind] = min(yyray(:,16,temp_num));%kper_min
        k_output(i,j,2) = min(yyray(ind,15,temp_num));%kparl
        k_output(i,j,3) = min(yyray(ind,18,temp_num));%ne
        k_output(i,j,4) = min(yyray(ind,1,temp_num));%r
        k_output(i,j,5) = min(yyray(ind,3,temp_num));%z
        k_output(i,j,7) = min(yyray(ind,4,temp_num));%z
        k_output(i,j,8) = min(yyray(ind,6,temp_num));%z
    end
%         [kper_min_down(i),ind] = min(yyray(:,16,i));
%         kper_down_kparl(i) = yyray(ind,15,i);
%         kper_down_r(i) = yyray(ind,1,i);
%         kper_down_z(i) = yyray(ind,3,i);
%         kper_down_B(i) = yyray(ind,17,i);
%         kper_down_ne(i) = yyray(ind,18,i);
end
%%
k_rr = k_output(:,:,4);
k_zz = k_output(:,:,5);
k_output(:,:,6) = interp2(rr',zz',rhorz',k_rr,k_zz);
%%
% f_ind = 1:2:num_f;
% theta_ind = (6:18) -1;
f_ind = 1:num_f;
theta_ind = 1:num_theta;

figure;

subplots(4,1,1);hold on;box on;
H1=plot(k_output(f_ind,theta_ind,6)',2*k_output(f_ind,theta_ind,1)'/100,'-o');
TitleLabels '' '\rho' 'k_{\perp}/cm^{-1}'
plot(k_output(f_ind,theta_ind,6),2*k_output(f_ind,theta_ind,1)/100,'-o');
legend(H1,string(f0(f_ind))+"GHz")

subplots(4,1,2);hold on;box on;
plot(k_output(f_ind,theta_ind,6)',-2*k_output(f_ind,theta_ind,2)'/100,'-o');
TitleLabels '' '\rho' 'k_{//}/cm^{-1}'
plot(k_output(f_ind,theta_ind,6),-2*k_output(f_ind,theta_ind,2)/100,'-o');

subplots(4,1,3);hold on;box on;
plot(k_output(f_ind,theta_ind,6)',-2*k_output(f_ind,theta_ind,7)'/100,'-o');
TitleLabels '' '\rho' 'k_{r}/cm^{-1}'
plot(k_output(f_ind,theta_ind,6),-2*k_output(f_ind,theta_ind,7)/100,'-o');

subplots(4,1,4);hold on;box on;
plot(k_output(f_ind,theta_ind,6)',-2*k_output(f_ind,theta_ind,8)'/100,'-o');
TitleLabels '' '\rho' 'k_{z}/cm^{-1}'
plot(k_output(f_ind,theta_ind,6),-2*k_output(f_ind,theta_ind,8)/100,'-o');


xlim([0.84 1.005]);
link x;