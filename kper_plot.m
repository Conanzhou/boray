
for i=1:size(yyray,3)
    kper_min(i) = min(yyray(:,16,i));
end
% figure;
plot(-1*(k_theta0-90),2*kper_min/100,'v');
