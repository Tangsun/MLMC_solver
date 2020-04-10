function main

clear all; close all; clc

%% system setup
M = 2;      %refinement factor
N_pilot = 1e5;    %pilot sample size
L_pilot = 8;     %initial levels with sample size

% Eps = [0.005, 0.01, 0.02, 0.05, 0.1];  %desired accuracies
Eps = [0.1 0.05 0.02 0.01];

%% Generate true value (estimation with large sample size)
x0 = 10; T = 1;
Q_true = x0*exp(-2*T + 0.5*T^2);

%% MLMC test results
test_result = mlmc_test(@linsys_lv, M, N_pilot, L_pilot);

figure(1); hold on
subplot(2, 2, 1); hold on
plot([1: L_pilot], log2(test_result.EQ(2: end)));
plot([1: L_pilot], log2(test_result.EP(2: end)));
xlabel('Level l'); ylabel('log_2 |mean|');
legend('Q_l', 'P_l = Q_l - Q_{l-1}');

subplot(2, 2, 2); hold on
plot([1: L_pilot], log2(test_result.varQ(2: end)));
plot([1: L_pilot], log2(test_result.varP(2: end)));
xlabel('Level l'); ylabel('log_2 var');
legend('Q_l', 'P_l = Q_l - Q_{l-1}');

subplot(2, 2, 3); hold on
plot([1: L_pilot], test_result.check(2: end), '--x');
xlabel('Level l'); ylabel('consistency check');

subplot(2, 2, 4); hold on
plot([0: L_pilot], test_result.kurt, '--x');
xlabel('Level l'); ylabel('kurtosis');

%% MLMC case
figure(2); subplot(2, 2, 1);
xlabel('Level l'); ylabel('N_l');
subplot(2, 2, 2);
xlabel('Mean Squared Error'); ylabel('Cost');
set(2,'DefaultAxesLineStyleOrder','--o|--x|--d');
mlmc_mse = zeros(length(Eps), 1); mlmc_cost = zeros(length(Eps), 1);
mlmc_var = zeros(length(Eps), 1); mlmc_err = zeros(length(Eps), 1);
mlmc_time = zeros(length(Eps), 1);
mc_mse = zeros(length(Eps), 1); mc_cost = zeros(length(Eps), 1);
mc_var = zeros(length(Eps), 1); mc_err = zeros(length(Eps), 1);
mc_time = zeros(length(Eps), 1);
for i = 1: length(Eps)
    [mlmc_sol, mc_sol] = mlmc(@linsys_lv, M, test_result, 3, Eps(i), Q_true, 500);
    subplot(2, 2, 1); hold on
    semilogy([0: mlmc_sol.L], mlmc_sol.Nl);
    mlmc_mse(i) = mlmc_sol.mse; mlmc_cost(i) = mlmc_sol.cost_pred;
    mlmc_var(i) = mlmc_sol.var; mlmc_err(i) = mlmc_sol.mu - Q_true;
    mlmc_time(i) = mlmc_sol.cost_act;
    mc_mse(i) = mc_sol.mse; mc_cost(i) = mc_sol.cost_pred;
    mc_var(i) = mc_sol.var; mc_err(i) = mc_sol.mu - Q_true;
    mc_time(i) = mc_sol.cost_act;
end
subplot(2, 2, 1);
legend('0.1', '0.05', '0.02', '0.01');
set(2,'DefaultAxesLineStyleOrder','-*|--*')
subplot(2, 2, 2); hold on
loglog(mlmc_mse, mlmc_cost);
loglog(mc_mse, mc_cost);
legend('MLMC', 'MC');
subplot(2, 2, 3); hold on
loglog(Eps, mlmc_var);
loglog(Eps, mc_var);
legend('MLMC', 'MC');
xlabel('\epsilon'); ylabel('variance');
subplot(2, 2, 4); hold on
loglog(Eps, abs(mlmc_err));
loglog(Eps, abs(mc_err));
legend('MLMC', 'MC');
xlabel('\epsilon'); ylabel('error');
figure(3);
subplot(2, 1, 1); hold on
loglog(Eps, mlmc_time, '-x');
loglog(Eps, mc_time, '-o');
legend('MLMC', 'MC');
xlabel('\epsilon'); ylabel('Time (N_{exp} = 500)');
subplot(2, 1, 2); hold on
loglog(Eps, mlmc_cost, '-x');
loglog(Eps, mc_cost, '-o');
xlabel('\epsilon'); ylabel('Pred Cost');
legend('MLMC', 'MC');
end