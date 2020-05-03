% Euler-Maruyama Scheme for SDEs
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
% W            Brownian motion
% Nsims        Number of simulation samples
% u            Control input
%
%%% Output
% t            Time
% out          State evolutions
%
%%% Notes
% "out" is of shape Nstates x Nsims x Nt

function [t, out] = euler_maruyama(X0, t0, dt, Nt, drift, sigma, W, Nsims)  
    Nstates = size(X0, 1); 
    out = zeros(Nstates, Nsims, Nt);
    t = zeros(Nt,1);
    t(1) = t0;
    %get Brownian motion increments
    dW = W(:, :, 2:end) - W(:, :, 1:end-1);
    %assign initial condition
    for i = 1:Nsims
        out(:, i, 1) = X0;
    end
    %simulate
    for i = 2:Nt
        %get drift
        out_drift = dt*drift(t(i-1), out(:, :, i-1)) + out(:, :, i-1);
        %get diffusion
        out_diffusion = dW(:, :, i-1) * sigma;
        %get simulation output
        out(:, :, i) = out_drift + out_diffusion';
        t(i) = t(i-1) + dt;
    end
end