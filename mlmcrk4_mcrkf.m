function [mlmc_sol, mc_sol] = mlmcrk4_mcrkf(objFcn, M, test_result, L_init, eps, Q_true, N_exp, reltol)
    %% Preprocessing
    %Read information at each level from test_result
    L = L_init; Nl = zeros(L+1, 1);
    L_update = 1;
    Vl = test_result.varP(1: L+1);
    El = test_result.EP(1: L+1);
    Cl = 2.^(test_result.gamma*(0: L))';
    Nl = ceil(2 * sqrt(Vl./Cl) * sum(sqrt(Vl.*Cl)) / eps^2);
    
    %Find the correct level number using linear regression parameter
    %Calculate theoretical prediction for Nl, Cl
    while L_update == 1
        range = -1:0;
        rem = max(El(L+1+range).*2.^(test_result.alpha*range)) / (2^test_result.alpha - 1);

        if rem > eps/sqrt(2)
            L       = L+1;
            Vl = test_result.varP(1: L+1);
            El = test_result.EP(1: L+1);
            Nl(L+1) = 0;

            Cl  = 2.^(test_result.gamma*(0:L))';
            Nl  = ceil(2 * sqrt(Vl./Cl) * sum(sqrt(Vl.*Cl)) / eps^2);
        else
            L_update = 0;
        end
    end
    
    %% MLMC estimation
    %Run MLMC experiments with L, Nl, Vl found earlier
    Y = zeros(N_exp, 1);
    tic;
    for i = 1: N_exp
        for l = 0: L
            [~, P_info] = feval(objFcn, l, Nl(l+1));
            Y(i) = Y(i) + P_info(1);
        end
    end
    mlmc_sol.cost_act = toc;
    mlmc_sol.cost_pred = 2*M*Nl(1);
    for i = 1: L
        mlmc_sol.cost_pred = mlmc_sol.cost_pred + Nl(i+1)*(2*M^i)*(1 + M);
    end
    mlmc_sol.L = L;
    mlmc_sol.Nl = Nl;
    
    mlmc_sol.mu = mean(Y);
    mlmc_sol.mse = sum((Y - Q_true).^2)/N_exp;
    mlmc_sol.var = var(Y);
    
    %% Equivalent standard MC estimation with rk4 at most coarse level
    %Run standard MC with equivalent computational cost
    Y_mc = zeros(N_exp, 1);
    mc_sol.cost_pred = 0;
    tic;
    for i = 1: N_exp
        X = zeros(100, 1);
        for j = 1: 100
            para = -2 + randn(1, 1);
            [~, X_i, ~, exec_i] = rkf45(@lin_sys, [0, 1], 10, 1e-3*reltol, reltol, -para);
            mc_sol.cost_pred = mc_sol.cost_pred + exec_i;
            X(j) = X_i(end);
        end
        Y_mc(i) = mean(X);
    end
    mc_sol.cost_act = toc;
    mc_sol.mu = mean(Y_mc);
    mc_sol.mse = sum((Y_mc - Q_true).^2)/N_exp;
    mc_sol.var = var(Y_mc);
end