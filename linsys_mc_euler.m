function Q_mc = linsys_mc_euler(l, N)
    
    T = 1; x0 = 10;
    M = 2; N0 = 2;
    
    nf = N0*M^(l+1);
    hf = T/nf;
    xf = x0*ones(1, N);
    
    a = -2 + randn(N, 1);
    A = sparse([1: N], [1: N], a, N, N);
    
    for j = 1: nf
        xf = xf + xf*A*hf;
    end
    Q_mc = mean(xf(end, :));
end