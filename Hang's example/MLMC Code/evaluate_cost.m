% sim_out is originally Nstates x Nsims x Nt
function cost = evaluate_cost(t, sim_out, Q)
cost = zeros(1, size(sim_out, 2));
% sim_out is now Nstates x Nt x Nsims
sim_out = permute(sim_out, [1 3 2]);
% Convert the angle so that is the one between pole and the vertical up
% position
sim_out(3, :, :) = wrapToPi(pi-sim_out(3, :, :));
for i = 1:size(sim_out, 3)
    sim_out_2d = sim_out(:, :, i);
    cost_temp = 0.5*sim_out_2d'*Q*sim_out_2d;
    cost_temp = diag(cost_temp)';
    cost(i) = trapz(t, cost_temp); 
end
end