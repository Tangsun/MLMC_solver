function [Q_info, P_info] = linsys_lv_euler(l, N)
    %Q_info, contains the mean of all objFcn samples, objFcn^2
    %P_info, contains the mean of all P_l, P_l^2, P_l^3, P_l^4
    Q_info = zeros(2, 1); P_info = zeros(4, 1);
    
    T = 1; x0 = 10;
    M = 2; N0 = 2;
    
    nf = N0*M^(l+1);
    hf = T/nf;
    xf = x0*ones(nf + 1, 1);
    
    if l > 0
        nc = nf/M;
        hc = T/nc;
        xc = x0*ones(nc + 1, 1);
    end
    
    for i = 1: N
        %generate random numbers for N samples
        a = -2 + randn(1, 1);
        
        for j = 1: nf
            xf(j+1) = xf(j) + a*xf(j)*hf;
        end
        Qf = xf(end);
        
        if l == 0
           Qc = 0;
        else
            for k = 1: nc
                xc(k+1) = xc(k) + a*xc(k)*hc;
            end
            Qc = xc(end);
        end
        Q_info(1) = Q_info(1) + Qf;
        Q_info(2) = Q_info(2) + Qf^2;
        P_info(1) = P_info(1) + Qf - Qc;
        P_info(2) = P_info(2) + (Qf - Qc)^2;
        P_info(3) = P_info(3) + (Qf - Qc)^3;
        P_info(4) = P_info(4) + (Qf - Qc)^4;
    end
    Q_info = Q_info/N;
    P_info = P_info/N;
end