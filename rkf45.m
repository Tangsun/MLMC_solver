function [T, X, err, exec] = rkf45(sys, tspan, x0, abstol, reltol, para)
    global a b c
    c = [0, 1/4, 3/8, 12/13, 1, 1/2];
    a = zeros(6, 5);
    a(2, 1) = 1/4;
    a(3, 1: 2) = [3/32, 9/32];
    a(4, 1: 3) = [1932/2197, -7200/2197, 7296/2197];
    a(5, 1: 4) = [439/216, -8, 3680/513, -845/4104];
    a(6, 1: 5) = [-8/27, 2, -3544/2565, 1859/4104, -11/40];
    b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;
        25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    
    if nargin == 3
        abstol = 1e-6;
    elseif nargin < 3
        printf('Not enough arguments');
    end
    
    exec = 0;
    h0 = 0.5; hmin = 3e-15;
    q = 1/2;
    T = tspan(1); X = x0;
    terminate = 0;
    x = X(:, 1); t = T(1);
    h = h0; err = [];
    while(~terminate)
        reltol_pass = 0;
        [x_rk5, e] = err_rk45(sys, h, t, x, para);
        if e < reltol*abs(x_rk5)
            reltol_pass = 1;
        elseif min(x_rk5) == 0
            reltol_pass = 1;
        end
        if max(e) <= abstol && reltol_pass == 1
            t = t + h;
            x = x_rk5;
            X = [X, x];
            T = [T, t];
            err = [err; e];
            h = 0.5;
        else
            s = (abstol/(2*max(e)))^0.25;
            h = s*h;
            if h < hmin
                h = hmin;
            elseif h > h0
                h = h0;
            end
        end
        if t >= tspan(2) - 1e-5
            terminate = 1;
        elseif t + h > tspan(2)
            h = tspan(2) - t;
        end
    end
end

function [x_rk5, e] = err_rk45(sys, h, t, x, para)
    global a b c

    K = zeros(size(x, 1), 6);
    for i = 1: 6
        if i == 1
            dx_T = zeros(1, size(x, 1));
        else
            dx_T = a(i, 1: (i-1))*K(:, 1: (i-1))';
        end
        K(:, i) = h*feval(sys, t + c(i)*h, x + dx_T', para);
    end
    x_rk5 = x + K*b(1, :)';
    x_rk4 = x + K*b(2, :)';
    e = abs(x_rk5 - x_rk4);
end