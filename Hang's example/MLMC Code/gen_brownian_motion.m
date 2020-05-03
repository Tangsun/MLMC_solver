%%% Generate Brownian motion
%
%%% last update: 28 October, 2019
%
%%% Inputs:
% Nsims        Number of simulation samples
% dt           Time step size
% Nt           Number of time steps
% Nstates      Number of states
%
%%% Output
% W            Brownian motion
%
%%% Notes
% "W" is of shape Nsims x Nstates x Nt

function W = gen_brownian_motion(Nsims, dt, Nt, Nstates)
    dW0 = zeros(Nsims, Nstates, Nt);
    dW0(:, :, 2:end) = sqrt(dt)*randn(Nsims, Nstates, Nt-1);
    W = cumsum(dW0, 3);
end