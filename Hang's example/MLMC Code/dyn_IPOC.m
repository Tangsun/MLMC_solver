function dXdt = dyn_IPOC(t, X, ufunc)
% Nonlinear model of an inverted pendulum on a cart
% Funtion inputs
% X     state vector (nstates x nsims)
% ufunc     control input
%
%%% Notes
% This function is for vectorized computation

% Free response
u = 0;

% Control
if nargin > 2
    u = ufunc(t, X); % u should ncontrols x nsims
end
    
% Setup
Nsims = size(X, 2);
dXdt = zeros(4, Nsims);

% Define parameters
M = 0.5;    %cart mass [kg]
m = 0.2;    %pendulum mass [kg]
b = 0.1;    %coefficient of friction for cart [N/m/sec]
I = 0.006;  %pendulum moment of inertia [kgm^2]
g = 9.8;    %gravitational acceleration [m/s^2]
l = 0.3;    %length to pendulum center of mass [m] (0.3 * 5)

% Scalars
p0 = m*l; 
p0g = p0 * g;
p02 = p0 * p0;
p1 = I+p0*l; % I + ml^2
p0gp1 = p0g / p1;

X(3, :) = pi + X(3, :);
% 1 x nsims
s3 = sin(X(3, :));
c3 = cos(X(3, :));
s3c3 = s3.*c3;
x42 = X(4, :).^2;
R = ((M+m)*p1 - p02 * c3.^2).^(-1);
Rx42 = R.*x42;
Rc3 = R.*c3;

% states
dXdt(1, :) = X(2, :);
dXdt(2, :) = p1*p0*Rx42.*s3 + p02*g*R.*s3c3 - p1*b*R.*X(2,:) + p1*R.*u;
dXdt(3, :) = X(4, :);
dXdt(4, :) = -p02*Rx42.*s3c3 - p02*p0gp1*Rc3.*s3c3 + p0*b*Rc3.*X(2,:) - p0*Rc3.*u - p0gp1*s3;
end