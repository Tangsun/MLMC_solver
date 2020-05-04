Nspl = mc_spl; 
eqL = [mlmc_spl{1}.L, mlmc_spl{2}.L, mlmc_spl{3}.L, mlmc_spl{4}.L, mlmc_spl{5}.L, mlmc_spl{6}.L];

N_exp = 500;
Y_mc = zeros(N_exp, 6);

for j = 1: 6
    for i = 1: N_exp
        Y_mc(i, j) = linsys_mc_euler(eqL(j), Nspl(j));
    end
end

bias = (mean(Y_mc) - Q_true).^2;
figure(5);
semilogy(Nspl, bias, 'x-');
xlabel('#Samples per experiement');
ylabel('bias^2');

figure(6); hold on
subplot(3, 2, 1);
histogram(Y_mc(:, 1)); xline(Q_true, 'r-');
xlabel('estimate Y'); ylabel('Occurance');
subplot(3, 2, 2);
histogram(Y_mc(:, 2)); xline(Q_true, 'r-');
xlabel('estimate Y'); ylabel('Occurance');
subplot(3, 2, 3);
histogram(Y_mc(:, 3)); xline(Q_true, 'r-');
xlabel('estimate Y'); ylabel('Occurance');
subplot(3, 2, 4);
histogram(Y_mc(:, 4)); xline(Q_true, 'r-');
xlabel('estimate Y'); ylabel('Occurance');
subplot(3, 2, 5);
histogram(Y_mc(:, 5)); xline(Q_true, 'r-');
xlabel('estimate Y'); ylabel('Occurance');
subplot(3, 2, 6);
histogram(Y_mc(:, 6)); xline(Q_true, 'r-');
xlabel('estimate Y'); ylabel('Occurance');