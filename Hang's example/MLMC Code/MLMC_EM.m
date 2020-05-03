% Execute MLMC algorithm with Euler-Maruyama
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
% np           Dimension of the objective
% 
%
%%% Output
% mean_est     MLMC mean estimate
%
%%% Notes
% "mean_est" is of shape Nt(1) x d

function P = MLMC_EM(X0, t0, dt, Nt, drift, sigma, N_samples, d, M, N_level, g, np)

% Variable set-up
P = zeros(np, N_level+1);

for i = 1:N_level+1
    if i == 1 % Level 0
        %get Brownian motion and use Euler-Maruyama to simulate
        [t, sim_out, ~] = MC_EM(X0, t0, dt(i), Nt(i), drift, sigma, N_samples(i), d);
        P(:, i) = mean(g(t, sim_out), 2);
    else  % Level 1 through N_level
        %fine evaluation
        [t, sim_out, W] = MC_EM(X0, t0, dt(i), Nt(i), drift, sigma, N_samples(i), d);
        %coarse evaluation
        W_coarse = brown_coarse_from_fine(W, M);
        [t_coarse, sim_out_coarse, ~] = MC_EM(X0, t0, dt(i-1), Nt(i-1), drift, sigma, N_samples(i), d, W_coarse);
        P(:, i) = mean(g(t, sim_out) - g(t_coarse, sim_out_coarse), 2);
    end
end

Cov = cov(P);
end