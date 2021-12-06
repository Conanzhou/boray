fkpar2 = @(kr)((kr * Br + kz * Bz + nphi / r * Bphi) / B).^2;
fkper2 = @(kr) (kr.^2 + kz^2 + nphi^2 / r^2) - fkpar2(kr);

% define the function to solve the initial kr
fDkr = @(kr)eps1 * (fkper2(kr) * c2 / w2).^2 - ...
    ((eps1 + eps3) * (eps1 - fkpar2(kr) * c2 / w2) - eps2^2) .* fkper2(kr) * c2 / w2 + ...
    eps3 * ((eps1 - fkpar2(kr) * c2 / w2).^2 - eps2^2);

krr = -1000:0.5:1000; krr = krr / 1;
subplot(121); plot(krr, fDkr(krr), krr, 0 .* krr, '--');
xlabel('k_r'); ylabel('D(\omega,k)');
subplot(122); plot(krr, sqrt(fkpar2(krr)), krr, sqrt(fkper2(krr)), '--');
xlabel('k_r'); ylabel('k'); legend('k_{||}', 'k_\perp'); legend('boxoff');
coeffa = eps1
eps2
eps3
fkpar2(1) * c2 / w2
coeffb = -((eps1 + eps3) * (eps1 - fkpar2(1) * c2 / w2) - eps2^2)
coeffc = eps3 * ((eps1 - fkpar2(1) * c2 / w2).^2 - eps2^2)
roots([coeffa, coeffb, coeffc])