function [Ql, Pl] = linsys_exp_euler(l, N)
    
    T = 1; x0 = 10;
    M = 2; N0 = 2;
    
    nf = N0*M^(l+1);
    hf = T/nf;
    xf = x0*ones(1, N);
    
    if l > 0
        nc = nf/M;
        hc = T/nc;
        xc = x0*ones(1, N);
    end
    
    a = -2 + randn(N, 1);
    A = sparse([1: N], [1: N], a, N, N);
    
    for j = 1: nf
        xf = xf + xf*A*hf;
    end
    Qf = xf(end, :);
    
    if l == 0
        Qc = 0;
    else
        for k = 1: nc
            xc = xc + xc*A*hc;
        end
        Qc = xc(end, :);
    end
    Pl = mean(Qf - Qc);
    Ql = mean(Qf);
end