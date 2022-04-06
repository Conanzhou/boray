%%
num_Freq = size(yyray,3)/2;

for i=1:num_Freq
    [kper_min_down(i),ind] = min(yyray(:,16,i));
    kper_down_kparl(i) = yyray(ind,15,i);
    kper_down_r(i) = yyray(ind,1,i);
    kper_down_z(i) = yyray(ind,3,i);
    kper_down_B(i) = yyray(ind,17,i);
    kper_down_ne(i) = yyray(ind,18,i);
end
for j=1:num_Freq
    i=j+num_Freq;
    [kper_min_up(j),ind] = min(yyray(:,16,i));
    kper_up_kparl(j) = yyray(ind,15,i);
    kper_up_r(j) = yyray(ind,1,i);
    kper_up_z(j) = yyray(ind,3,i);
    kper_up_B(j) = yyray(ind,17,i);
    kper_up_ne(j) = yyray(ind,18,i);
end
%%
figure;
subplots(4,1,1);
plot(f0(1:num_Freq),2*kper_min_down /100,'-o');hold on;
plot(f0(1:num_Freq),2*kper_min_up /100,'-s');
TitleLabels '' 'Freq/GHz' 'k_{\perp}/cm^{-1}'
subplots(4,1,2);
plot(f0(1:num_Freq),2*kper_down_kparl /100,'-o');hold on;
plot(f0(1:num_Freq),2*kper_up_kparl /100,'-s');
TitleLabels '' 'Freq/GHz' 'k_{//}/cm^{-1}'
subplots(4,1,3);
plot(f0(1:num_Freq),kper_down_B,'-o');hold on;
plot(f0(1:num_Freq),kper_up_B,'-s');
TitleLabels '' 'Freq/GHz' 'B/T'
subplots(4,1,4);
plot(f0(1:num_Freq),kper_down_ne/1e19,'-o');hold on;
plot(f0(1:num_Freq),kper_up_ne/1e19,'-s');
TitleLabels '' 'Freq/GHz' 'ne/\times 19 m^{-3}'

xlim([51 66]);
link x;
