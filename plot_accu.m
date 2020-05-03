eps = Eps;
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