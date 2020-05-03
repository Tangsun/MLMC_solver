function Nl_opt = Optimal_Allocation(V, CL, M, Var_target, N_level)
% Calculate the optimal sample allocation for MLMC
C_hat = zeros(N_level+1, 1);
Nl_opt = zeros(N_level+1, 1);
lambda = 0;

for i = 1:N_level+1
    C_hat(i) = CL/M^(N_level+1-i);
    lambda = lambda + sqrt(V(i)*C_hat(i));
end

lambda = lambda*Var_target^(-2);

for i = 1:N_level+1
    Nl_opt(i) = sqrt(V(i)/C_hat(i));
end

Nl_opt = lambda*Nl_opt;
end