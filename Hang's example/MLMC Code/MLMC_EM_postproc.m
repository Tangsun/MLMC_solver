% Execute MLMC algorithm with Euler-Maruyama with post-processing
%
%%% last update: 29 October, 2019
%
%%% Inputs:
% X0           Initial condition
% dt           A (N_level+1) by 1 array with each element being the time
%              step size for the corresponding level
% Nt           A (N_level+1) by 1 array with each element being the number
%              of time steps for the corresponding level
% drift        Drift
% sigma        Diffusion
% N_samples    A (N_level+1) by 1 arrary containing the sample allocation
% d            Number of states
% M            Refinement factor between two neighboring levels
% N_level      Number of MLMC level
% g            Post-processing function to evaluate value/cost
% np           Dimension of the output of g
% 
%
%%% Output
% mean_est     MLMC mean estimate
%
%%% Notes
% "mean_est" is of shape Nt(1) x d

function [P, P_var, P_level_fine, P_var_level_fine, g_out, g_out_coarse] = MLMC_EM_postproc(X0, t0, dt, Nt, drift, sigma, N_samples, d, M, N_level, g, np)

% Variable set-up
P = zeros(np, N_level+1);
P_var = zeros(np, N_level+1);
P_level_fine = zeros(np, N_level);
P_var_level_fine = zeros(np, N_level);
g_out = cell(N_level+1);
g_out_coarse = cell(N_level, 1);

for i = 1:N_level+1
    if i == 1 % Level 0
        %get Brownian motion and use Euler-Maruyama to simulate
        [t, sim_out, ~] = MC_EM(X0, t0, dt(i), Nt(i), drift, sigma, N_samples(i), d);
        P(:, i) = mean(g(t, sim_out), 2);
        %post-processing calculations
        P_var(:, i) = var(g(t, sim_out), 0, 2);
        g_out{i} = sim_out;
    else  % Level 1 through N_level
        %fine evaluation
        [t, sim_out, W] = MC_EM(X0, t0, dt(i), Nt(i), drift, sigma, N_samples(i), d);
        %post-processing calculations
        P_level_fine(:, i-1) = mean(g(t, sim_out), 2);
        P_var_level_fine(:, i-1) = var(g(t, sim_out), 0, 2);
        g_out{i} = sim_out;
        %coarse evaluation
        W_coarse = brown_coarse_from_fine(W, M);
        [t_coarse, sim_out_coarse, ~] = MC_EM(X0, t0, dt(i-1), Nt(i-1), drift, sigma, N_samples(i), d, W_coarse);
        P(:, i) = mean(g(t, sim_out) - g(t_coarse, sim_out_coarse), 2);
        %post_processing calculations
        P_var(:, i) = var(g(t, sim_out) - g(t_coarse, sim_out_coarse), 0, 2);
        g_out_coarse{i-1} = sim_out_coarse;
    end
end
end