%%% Get Coarse Brownian Mition from Fine Brownian Motion
%
%%% last update: 28 October, 2019
%
%%% Inputs:
% Nsims        Fine Brownian motion
% M            Refinement factor between two neighboring levels
%
%%% Output
% coarse       Corase Brownian motion

function coarse = brown_coarse_from_fine(fine, M)
    coarse = fine(:, :, 1:M:size(fine,3));
end