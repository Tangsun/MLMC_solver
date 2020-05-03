function Ql = linsys_mc_rk4(l, N)

    T = 1; x0 = 10;
    M = 2; N0 = 2;
    
    nf = N0*M^(l+1);
    hf = T/nf;
    xf = x0*ones(1, N);
    
    %generate random numbers for N samples
    a = -2 + randn(1, 1);
    A = sparse([1: N], [1: N], a, N, N);
    
    for j = 1: nf
        k1 = hf*xf*A;
        k2 = hf*(xf + k1/2)*A;
        k3 = hf*(xf + k2/2)*A;
        k4 = hf*(xf + k3)*A;
        xf = xf + 1/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    Ql = mean(xf(end, :));
end