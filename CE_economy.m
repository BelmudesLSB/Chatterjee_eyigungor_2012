%%% MATLAB:
% This code Solves the Chatterje and Eyigungor model using open mp.
 
clc;
clear;

%% Step 0 : we need to compile the MEX code.

mex -g main.cpp economy.cpp ggq_algorithm.cpp econ_functions.cpp math_functions.cpp aux_functions.cpp
%mex CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" main.cpp ArellanoEconomy.cpp Utils.cpp

%% Step 1: Assuming you have compiled the C++ code into a MEX file named "main_mex"

% Set up the parameters for the economy:
tic;
params.beta = 0.95402;       % discount factors
params.gamma = 2;           % risk aversion

% Income process parameters:
params.rho = 0.948503;      % persistence of income
params.sigma_e = 0.027092;  % standard deviation of income shocks
params.sigma_m = 0.003;     % standard deviation of income shocks
params.t = 3;               % number of standard deviations for Tauchen (1986)
params.d_0 = -0.18819;      % Output loss at default parameter
params.d_1 = 0.24558;       % Output loss at default parameter

% Grid Parameters:
params.ny = 21;             % Y grid points
params.nb = 152;            % B grid points
params.b_min = -0.8;        % Minimum bond holdings
params.b_max = 0.0;         % Maximum bond holdings
params.m_bar = 0.006;       % Maximum value for m

% Convergence parameters:
params.max_iter = 50000;    % Maximum number of iterations
params.tol = 1e-6;          % Tolerance for convergence

% Term structure parameters:
params.r = 0.01;            % Risk free interest rate
params.lambda = 0.05;       % Reciprocal of average maturity
params.z = 0.03;            % Coupon payments
params.xi = 0.035;          % Probability of reentry
params.eta_q = 0.995;       % weight on the old bond price
params.eta_w = 0.9;         % weight on the new continuation value
params.eta_vd = 0;          % weight on the new value of default


% Call the MEX function to solve the model
output = main(params);
toc;
%%

% Retrieve the output arrays
ygrid = output.Ygrid;
bgrid = output.Bgrid;
p = output.P;
v = output.V;
v_r = output.V_r;
v_d = output.V_d;
q = output.Q;
b_p = output.B_p;
d_p = output.D_p;
toc;

% Transform everything from vectors into matrices (n_y x n_b):
Y = ygrid;
B = bgrid;
P = (reshape(p, [params.ny, params.ny]))';
V = (reshape(v, [params.nb, params.ny]))';
V_r = (reshape(v_r, [params.nb, params.ny]))';
Q = (reshape(q, [params.nb, params.ny]))';
D_p = (reshape(d_p, [params.nb, params.ny]))';
B_p = (reshape(b_p, [params.nb, params.ny]))';
ny = params.ny;
nb = params.nb;