%%
% num_Freq = size(yyray,3)/2;

% k_output = zeros(num_f,num_theta,6);
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
% k_output(:,:,6) = sqrt((k_output(:,:,1)-R0).^2+(k_output(:,:,3)-Z0).^2)/a;%rho
%%
save([savepath,'k_output_DBS_H_up.mat'],'k_output');
%%
figure;hold on;box on;
for j = 1:num_theta
    plot(f0,2*k_output(:,j,1)/100,'-o');
end
legend(string(90-k_theta0(1:num_theta))+"^o");
TitleLabels '' 'Freq/GHz' 'k_{\perp}/cm^{-1}'
xlim([49 70]);
saveas(gcf,[savepath,'kper~f.fig']);
figure;hold on;box on;
for j = 1:num_theta
    plot(f0,2*k_output(:,j,2)/100,'-o');
end
legend(string(90-k_theta0(1:num_theta))+"^o");
TitleLabels '' 'Freq/GHz' 'k_{//}/cm^{-1}'
xlim([49 70]);
saveas(gcf,[savepath,'kparl~f.fig']);
%%
figure;hold on;box on;
for j = 1:num_f
    plot(90-k_theta0,2*k_output(j,:,1)/100,'-o');
end
legend(string(f0(1:num_f))+"GHz");
TitleLabels '' 'theta' 'k_{\perp}/cm^{-1}'
saveas(gcf,[savepath,'kper~theta.fig']);
figure;hold on;box on;
for j = 1:num_f
    plot(90-k_theta0,2*k_output(j,:,2)/100,'-o');
end
legend(string(f0(1:num_f))+"GHz");
TitleLabels '' 'theta' 'k_{//}/cm^{-1}'


%%
k_rr = k_output(:,:,4);
k_zz = k_output(:,:,5);
k_output(:,:,6) = interp2(rr',zz',rhorz',k_rr,k_zz);
figure;hold on;box on;
for j = 1:num_theta
    plot(k_output(:,j,6),2*k_output(:,j,1)/100,'-o');
end
legend(string(90-k_theta0(1:num_theta))+"^o");
TitleLabels '' '\rho' 'k_{\perp}/cm^{-1}'
for j = 1:num_f
    plot(k_output(j,:,6),2*k_output(j,:,1)/100,'-o');
end

%%
figure;axes;hold on
% H=plot(k_output(:,:,6)',2*k_output(:,:,1)'/100,'-o');

plot(k_output(:,1:12,6),2*k_output(:,1:12,1)/100,'-');
plot(k_output(:,23:-1:12,6),2*k_output(:,23:-1:12,1)/100,'-');
H=plot(k_output(:,1:12,6)',2*k_output(:,1:12,1)'/100,'-o');
plot(k_output(:,23:-1:12,6)',2*k_output(:,23:-1:12,1)'/100,'-*');
plot([0.9 0.9],[0 25]);
legend(H,string(f0)+"GHz")
TitleLabels 'ne=4.5e19,Bt=1.8T,Up Window' '\rho' 'k_{\perp}/cm^{-1}'
xlim([0.2 1.01]);

% plot(k_output(:,1:12,6),2*k_output(:,1:12,1)/100,'-');
% plot(k_output(:,20:-1:12,6),2*k_output(:,20:-1:12,1)/100,'-');
% H=plot(k_output(:,1:12,6)',2*k_output(:,1:12,1)'/100,'-o');
% plot(k_output(:,20:-1:12,6)',2*k_output(:,20:-1:12,1)'/100,'-');
% plot([0.9 0.9],[0 25]);
% legend(H,string(f0)+"GHz")
% TitleLabels 'ne=4.5e19,Bt=1.8T,Middle Window' '\rho' 'k_{\perp}/cm^{-1}'
% xlim([0.8 1.01]);
%%
figure;axes;hold on
% plot(k_output(:,1:11,6),2*k_output(:,1:11,1)/100,'-o');

% plot(k_output(:,21:-1:11,6),2*k_output(:,21:-1:11,2)/100,'-o');
% plot(k_output(:,1:11,6)',2*k_output(:,1:11,2)'/100,'-o');
% H=plot(k_output(:,21:-1:11,6)',2*k_output(:,21:-1:11,2)'/100,'-o');
% plot([0.9 0.9],[-6 6]);
% legend(H,string(f0)+"GHz")

plot(k_output(:,1:12,6),2*k_output(:,1:12,2)/100,'-');
plot(k_output(:,23:-1:12,6),2*k_output(:,23:-1:12,2)/100,'-');
H=plot(k_output(:,1:12,6)',2*k_output(:,1:12,2)'/100,'-o');
plot(k_output(:,23:-1:12,6)',2*k_output(:,23:-1:12,2)'/100,'-');
plot([0.9 0.9],[-6 6]);
legend(H,string(f0)+"GHz")
TitleLabels 'ne=4.5e19,Bt=1.8T,Up Window' '\rho' 'k_{\perp}/cm^{-1}'
xlim([0.2 1.01]);

% plot(k_output(:,1:12,6),2*k_output(:,1:12,2)/100,'-');
% plot(k_output(:,20:-1:12,6),2*k_output(:,20:-1:12,2)/100,'-');
% H=plot(k_output(:,1:12,6)',2*k_output(:,1:12,2)'/100,'-o');
% plot(k_output(:,20:-1:12,6)',2*k_output(:,20:-1:12,2)'/100,'-');
% plot([0.9 0.9],[-6 6]);
% legend(H,string(f0)+"GHz")
% TitleLabels 'ne=4.5e19,Bt=1.8T,Middle Window' '\rho' 'k_{//}/cm^{-1}'
% xlim([0.8 1.01]);