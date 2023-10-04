%% Step 0: Load the results and the parameters from the CE_economy file;
clear all
clc
load('C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\4\MATLAB\Parameters.mat')
load('C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\4\MATLAB\Results_calibrated_model.mat')


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
q_10 = 0;
W_10 = 0;

%% Start the algorithm

% Check for feasible consumptions
for x=1:params.nb
    if c(x)+params.m_bar>0
        X_tilde(end+1) = x;
    end
end
 
Exit = false;

% If X_tilda is empty.
if isempty(X_tilde)
    Exit = true;
    m_grid = [-params.m_bar, params.m_bar];
    d_grid = 1;
    b_grid = nan;
    q_10 = 0;
    W_10 = output.Vd(i);
end

% If X_tilda is non-empty: Get the optimal policy associated with m_bar.
if (~Exit)
    temp1 = -inf;
    x_1 = nan;
    for x=1:params.nb
        temp2 = utility(c(x)+params.m_bar, params.gamma)+W(x);
        if temp2 > temp1
            temp1 = temp2;
            x_1 = x;
        elseif temp2 == temp1
            if c(x) < c(x_1)
                temp1 = temp2;
                x_1 = x;
            end
        end
    end
end

if (~Exit)
    x_i = x_1;
    m_i = params.m_bar;
    % Create the objects to fill:
    m_grid = [];
    d_grid = [];
    b_grid = [];
end

% If X_tilda is not empty.
while ~Exit
    % First check if the policy under scrutiny is not preferred to default.
    if (utility(c(x_i)+params.m_bar, params.gamma)+W(x_i)<=output.Vd(i))
        d_grid = [1, d_grid];
        m_grid = [-params.m_bar, m_i, m_grid];
        b_grid = [nan, b_grid];
        q_10 = 0;
        W_10 = Vd(i)*cdf_truncated(m_i, -params.m_bar, params.m_bar, params.sigma_m);
        Exit = true;
        break;
    end
    % Default is not prefered.
    X_tilde = []
    for x=1:params.nb
        if c(x)>c(x_1)
            X_tilde(end+1) = x;
        end
    end
    if isempty(X_tilde)
        % There is no other policy to dispute the current one:
        % FUNCTION_LAST_UPDATE
        if (c(x_i)-params.m_bar)>0
            if utility(c(x_i)-params.m_bar, params.gamma)+W(x_1)>=output.Vd(i)
                m_grid = [-params.m_bar, m_i];
                d_grid = [0, d_grid];
                b_grid = [x_i, b_grid];
                W_10 = W_10 + (quadgk(@(m)(utility(c(x_i)+m, params.gamma)+W(x_1).*normpdf(m,0,params.sigma_m)),-m_bar, m_i))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m));
                q_10 = q_10 + (1/(1+params.r)) * (params.lambda + (1-params.lambda)*(z+ouput.q(i*nb+x_i))) * ((normcdf(m_i, 0, params.sigma_m)-normcdf(-m_bar, 0, params.sigma_m ))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m)));
                Exit = true;
                break;
            else
                m_indf = ((output.Vd(i) - W(x_i))*(1 - params.gamma))^(1/(1-params.gamma)) - c(x_i);
                m_grid = [-params.m_bar, m_inf, m_i, m_grid];
                b_grid = [nan, x_i, b_grid];
                d_grid = [1, 0, d_grid];
                % Not default.
                W_10 = W_10 + (quadgk(@(m)(utility(c(x_i)+m, params.gamma)+W(x_1).*normpdf(m,0,params.sigma_m)),m_indf, m_i))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m));
                q_10 = q_10 + (1/(1+params.r)) * (params.lambda + (1-params.lambda)*(z+ouput.q(i*nb+x_i))) * ((normcdf(m_i, 0, params.sigma_m)-normcdf(m_indf, 0, params.sigma_m))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m)));
                % Default.
                W_10 = W_10 + cdf_truncated(m_indf, -params.m_bar, params.m_bar, params.sigma_m) * output.Vd(i);
            end
        else
            m_indf = ((output.Vd(i) - W(x_i))*(1 - params.gamma))^(1/(1-params.gamma)) - c(x_i);
            m_grid = [-params.m_bar, m_inf, m_i, m_grid];
            b_grid = [nan, x_i, b_grid];
            d_grid = [1, 0, d_grid];
            % Not default.
            W_10 = W_10 + (quadgk(@(m)(utility(c(x_i)+m, params.gamma)+W(x_1).*normpdf(m,0,params.sigma_m)),m_indf, m_i))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m));
            q_10 = q_10 + (1/(1+params.r)) * (params.lambda + (1-params.lambda)*(z+ouput.q(i*nb+x_i))) * ((normcdf(m_i, 0, params.sigma_m)-normcdf(m_indf, 0, params.sigma_m))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m)));
            % Default.
            W_10 = W_10 + cdf_truncated(m_indf, -params.m_bar, params.m_bar, params.sigma_m) * output.Vd(i);
            Exit = True;
            break
        end
        
    else
        m_hat = -inf;
        x_hat = nan;
        for x=1:params.nb
            if c(x)>c(x_1)
                m_candidate = Find_m_crossing(params.gamma, W(x_i), W(x), c(x_i), c(x), m_i);
                if m_candidate > aux
                    m_hat = m_candidate;
                    x_hat = x;
                elseif m_candidate == aux
                    if c(x)>c(x_hat)
                        m_hat = m_candidate;
                        x_hat = x;
                    end
                end
            end        
        end

        if m_hat>-params.m_bar
            if utility(c(x_i) + m_hat, params.gamma) + W(x_i) >= output.Vd(i)
                m_grid = [m_hat, m_i];
                d_grid = [0, d_grid];
                b_grid = [x_i, b_grid];
                W_10 = W_10 + (quadgk(@(m)((utility(c(x_i)+m, params.gamma)+W(x_1))*normpdf(m,0,params.sigma_m)),m_hat, m_i))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m));
                q_10 = q_10 + (1/(1+params.r)) * (params.lambda + (1-params.lambda)*(z+ouput.q(i*nb+x_i))) * (normcdf(m_i, 0, params.sigma_m - normcdf(m_hat, 0, params.sigma_m))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m)));
                x_i = x_hat;
                m_i = m_hat;
            else
                m_indf = ((output.Vd(i) - W(x_i))*(1 - params.gamma))^(1/(1-params.gamma)) - c(x_i);
                m_grid = [-params.m_bar, m_inf, m_i, m_grid];
                b_grid = [nan, x_i, b_grid];
                d_grid = [1, 0, d_grid];
                % Not default.
                W_10 = W_10 + (quadgk(@(m)((utility(c(x_i)+m, params.gamma)+W(x_1))*normpdf(m,0,params.sigma_m)),m_indf, m_i))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m));
                q_10 = q_10 + (1/(1+params.r)) * (params.lambda + (1-params.lambda)*(z+ouput.q(i*nb+x_i))) * ((normcdf(m_i, 0, params.sigma_m)-normcdf(m_indf, 0, params.sigma_m))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m)));
                % Default.
                W_10 = W_10 + cdf_truncated(m_indf, -params.m_bar, params.m_bar, params.sigma_m) * output.Vd(i);
                Exit = True;
                break 
            end
        else %m_hat<=-params.m_bar
            if utility(c(x_i) -params.m_bar, params.gamma) + W(x_i) >= output.Vd(i)
                m_grid = [-params.m_bar, m_i];
                d_grid = [0, d_grid];
                b_grid = [x_i, b_grid];
                W_10 = W_10 + (quadgk(@(m)((utility(c(x_i)+m, params.gamma)+W(x_1))*normpdf(m,0,params.sigma_m)),-m_bar, m_i))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m));
                q_10 = q_10 + (1/(1+params.r)) * (params.lambda + (1-params.lambda)*(z+ouput.q(i*nb+x_i))) * (normcdf(m_i, 0, params.sigma_m - normcdf(-m_bar, 0, params.sigma_m ))/(normcdf(params.m_bar, 0, params.sigma_m)-normcdf(-params.m_bar, 0, params.sigma_m)));
                Exit = true;
                break;
            else 
            end
        end
    end  
end










%% Functions

% Output loss:
function phi = phi_f(i, ygrid, d0, d1)
     phi = max(0, d0*ygrid(i)+d1*ygrid(i)^2);
end

% Utility
function u = utility(c, gamma)
    u = (1-gamma)*c^(1-gamma);
end

% Cumulative prob truncated
function cum_prob = cdf_truncated(x, a, b, sigma)
    den = (normcdf(b, 0, sigma) - normcdf(a, 0, sigma));
    if x>=a && x<=b 
        cum_prob = normcdf(x, 0, sigma)/den;
    else
        cum_prob = nan;
    end
end

% Finding the m that satisfies the Step 2.
function m_star = Find_m_crossing(gamma, W_i, W_hat, c_i, c_hat, m_i)
    objective = @(m) utility(c_i + m, gamma) + W_i - utility(c_hat + m, gamma) + W_hat;
    m_star = fsolve(objective, m_i);
end


