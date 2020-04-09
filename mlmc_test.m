function test_result = mlmc_test(objFcn, M, N, L_pilot)
    rng('default');     % reset random number generator
    
    EQlv    = zeros(L_pilot + 1, 1);
    varQlv  = zeros(L_pilot + 1, 1);
    EPlv    = zeros(L_pilot + 1, 1);
    varPlv  = zeros(L_pilot + 1, 1);
    kurlv   = zeros(L_pilot + 1, 1);
    checklv = zeros(L_pilot + 1, 1);
    cost    = zeros(L_pilot + 1, 1);
    
    for l = 0: L_pilot
        fprintf('Pilot sample calculated at Level: l = %d\n', l);
        tic;
        [Q_info, P_info] = feval(objFcn, l, N);
        cost(l+1) = toc;
        EPlv(l+1)   = P_info(1);
        EQlv(l+1)   = Q_info(1);
        varPlv(l+1) = P_info(2) - P_info(1)^2;
        varQlv(l+1) = max(Q_info(2) - Q_info(1)^2, 1e-12);
        kurlv(l+1)  = (P_info(4) - 4*P_info(3)*P_info(1) + 6*P_info(2)*P_info(1)^2 - ...
            3*P_info(1)^4)/(P_info(2) - P_info(1)^2)^2;
        if l == 0
            checklv(l+1) = 0;
        else
            checklv(l+1) = abs(EPlv(l+1) + EQlv(l) - EQlv(l+1))/( 3*(sqrt(varPlv(l+1))...
                + sqrt(varQlv(l)) + sqrt(varQlv(l+1)) )/sqrt(N) );
        end
    end
    
    % estimate alpha, beta, gamma
    range = max(2, floor(0.4*(L_pilot+1))): (L_pilot+1);
    fprintf('\nestimate MLMC Theorem parameter based on linear regression:\n');
    pa = polyfit((range - 1)', log2(abs(EPlv(range))), 1);
    alpha = -pa(1); c1 = pa(2);
    fprintf('alpha = %f, c_1 = %f, for weak convergence\n', alpha, c1);
    pb = polyfit((range - 1)', log2(abs(varPlv(range))), 1);
    beta = -pb(1); c2 = pb(2);
    fprintf('beta = %f, c_2 = %f, for variance decay\n', beta, c2);
    pc = polyfit((range - 1)', log2(abs(cost(range))), 1);
    gamma = pc(1); c3 = cost(end)/(2^(gamma*(L_pilot)));
    fprintf('gamma = %f, c_3 = %f, for cost\n', gamma, c3);
    
    if max(checklv) > 1
        fprintf('\nWARNING: maximum consistency error = %f \n', max(check1v));
        fprintf('indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n\n');
    end
    
    if kurlv(end) > 100
        fprintf('WARNING: kurtosis on finest level = %f \n', kurlv(end));
        fprintf('indicates MLMC correction dominated by a few rare paths; \n');
    end
    
    test_result.EQ     = EQlv;
    test_result.varQ   = varQlv;
    test_result.EP     = EPlv;
    test_result.varP   = varPlv;
    test_result.kurt   = kurlv;
    test_result.check  = checklv;
    test_result.cost   = cost;
    test_result.alpha = alpha;
    test_result.beta = beta;
    test_result.gamma = gamma;
end