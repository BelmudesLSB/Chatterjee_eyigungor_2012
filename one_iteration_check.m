%% Step 0: Load the results and the parameters from the CE_economy file;

load('Results_calibrated_model.mat')
load('Parameters.mat');

%% Step 1: 

% Fix (i,b):

i = 10;
b = 100;

% Given price define the vector for consumption:

aux_1 = zeros(params.nb);
aux_2 = zeros(params.nb);
for x=1:params.nb
    aux_1(x) = -output.Q(i * params.nb + x) * (output.bgrid(x) - (1-params.lambda) * output.bgrid(b)) + (params.lambda + params.z * (1-params.lambda)) * output.bgrid(b);
    aux_2(x) = ouput.W(i * params.nb + x);
end

c = aux;
W = aux;

X_tilde = [];




