%% Step 0: Load the results and the parameters from the CE_economy file;
clear all
clc
cd 'C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\MATLAB'
load('C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\MATLAB\Parameters.mat')
load('C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\MATLAB\Results_calibrated_model.mat')


%% Step 1: 

% Fix (i,b):
i = 10;
b = 100;

% Given price define the vector for consumption:

aux_1 = zeros(params.nb);
aux_2 = zeros(params.nb);

for x=1:params.nb
    aux_1(x) = -output.Q(i * params.nb + x) * (output.bgrid(x) - (1-params.lambda) * output.bgrid(b)) + (params.lambda + params.z * (1-params.lambda)) * output.bgrid(b);
    aux_2(x) = output.W(i * params.nb + x);
end

c = aux_1;
W = aux_2;
X_tilde = [];

phi(2)
f1(3)

%%
for x=1:params.nb
    if c(x)>0
        X_tilde(end+1) = x;
    end
end

Exit = false;

if isempty(X_tilde)
    Exit = true;
    m_grid = [-params.m_bar, params.m_bar];
    d_grid = 1;
    b_grid = nan;
    c_grid = output.ygrid(i) - phi(i) - params.m_bar;
end

while ~Exit


end



