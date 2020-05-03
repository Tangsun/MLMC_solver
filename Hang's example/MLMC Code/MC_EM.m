% Monte Carlo sampling of a stochastic differential equation with Euler-Maruyama
%
%%% last update: 29 October, 2019
%
%%% Inputs:
% X0           Initial condition
% t0           Initial time
% dt           Time step size
% Nt           Number of time steps
% drift        Drift
% sigma        Diffusion
% N_samples    Sample size
% d            Number of states
% W            Brownian motion (if provided)
%
%%% Output
% t_out        Time
% sim_out      Simulation outputs     
% W_out        Brownian motion
%
%%% Notes
% "sim_out" is of shape Nstates x Nsims x Nt

function [t_out, sim_out, W_out] = MC_EM(X0, t0, dt, Nt, drift, sigma, N_samples, d, W)
% Get Brownian motion
if nargin == 8
    %generate Brownian motion for all states for all samples
    W_out = gen_brownian_motion(N_samples, dt, Nt, d);
else
    %use the provided Brownian motion
    W_out = W;
end

% Simulate with Euler-Maruyama
[t_out, sim_out] = euler_maruyama(X0, t0, dt, Nt, drift, sigma, W_out, N_samples);

end