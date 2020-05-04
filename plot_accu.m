eps = Eps;
%% plot results
figure(2);
for i = 1: length(Eps)
    subplot(2, 2, 1); hold on
    semilogy([0: mlmc_spl{i}.L], mlmc_spl{i}.spl);
end
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('0.03', '0.025', '0.02', '0.015', '0.01', '0.005');
xlabel('Level l'); ylabel('N_l');
subplot(2, 2, 2); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Mean Squared Error'); ylabel('Cost');
plot(mlmc_mse, mlmc_cost, '-x');
plot(mc_mse, mc_cost, '-o');
legend('MLMC', 'MC');
subplot(2, 2, 3); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
loglog(Eps, mlmc_time, '-x');
loglog(Eps, mc_time, '-o');
legend('MLMC', 'MC');
xlabel('\epsilon'); ylabel('Time (N_{exp} = 500)');
subplot(2, 2, 4); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
loglog(Eps, mlmc_cost, '-x');
loglog(Eps, mc_cost, '-o');
xlabel('\epsilon'); ylabel('Pred Cost');
legend('MLMC', 'MC');
figure(3);
subplot(2, 2, 1); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
plot(Eps, mlmc_var, '-x');
plot(Eps, mc_var, '-o');
legend('MLMC', 'MC');
xlabel('\epsilon'); ylabel('variance');
subplot(2, 2, 2); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
plot(Eps, mlmc_err, '-x');
plot(Eps, mc_err, '-o');
legend('MLMC', 'MC');
xlabel('\epsilon'); ylabel('bias^2');
subplot(2, 2, 3); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
plot(mlmc_time, mlmc_var, '-x');
plot(mc_time, mc_var, '-o');
plot(mlmc_cost, mlmc_var, '--x');
plot(mc_cost, mc_var, '--o');
legend('MLMC act', 'MC act', 'MLMC pred', 'MC pred');
xlabel('Act cost'); ylabel('variance');
subplot(2, 2, 4); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
plot(mlmc_time, mlmc_err, '-x');
plot(mc_time, mc_err, '-o');
plot(mlmc_cost, mlmc_err, '--x');
plot(mc_cost, mc_err, '--o');
legend('MLMC act', 'MC act', 'MLMC pred', 'MC pred');
xlabel('Act cost'); ylabel('bias^2');

figure(4);
subplot(1, 2, 1); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
plot(mlmc_cost, mlmc_var + mlmc_err.^2, '-x');
plot(mlmc_cost, mlmc_mse, '-.');
plot(mlmc_cost, eps.^2, '--x');
plot(mc_cost, mc_var + mc_err.^2, '-o');
plot(mc_cost, mc_mse, '-*');
plot(mc_cost, eps.^2, '--o');
legend('MLMC var+bias', 'MLMC mse', 'MLMC \epsilon^2', 'MC var+bias', 'MC mse', 'MC \epsilon^2');
xlabel('Pred cost'); ylabel('accuracy');
subplot(1, 2, 2); hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
plot(mlmc_time, mlmc_var + mlmc_err.^2, '-x');
plot(mlmc_time, mlmc_mse, '-.');
plot(mlmc_time, eps.^2, '--x');
plot(mc_time, mc_var + mc_err.^2, '-o');
plot(mc_time, mc_mse, '-*');
plot(mc_time, eps.^2, '--o');
legend('MLMC var+bias', 'MLMC mse', 'MLMC \epsilon^2', 'MC var+bias', 'MC mse', 'MC \epsilon^2');
xlabel('Act cost'); ylabel('accuracy');