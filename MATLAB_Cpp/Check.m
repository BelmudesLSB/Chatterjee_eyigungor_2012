%% Step: Load the results and the parameters from the CE_economy file;
clear all, clc
load('C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\4\MATLAB\Parameters.mat')
load('C:\Users\belmu\OneDrive\Escritorio\Repositories\Check\4\MATLAB\Results_calibrated_model.mat')

%% Step: Unload parameters:

gamma = params.gamma;
beta = params.beta;
m_bar = params.m_bar;
sigma_m = params.sigma_m;
d_0 = params.d_0;
d_1 = params.d_1;
nb = params.nb;
ny = params.ny;
lambda = params.lambda;
z = params.z;
r = params.r;

P = reshape(output.pgrid, ny, ny)';

Q_0 = reshape(output.Q, nb, ny)'; % Q(i,b) the price of debt b issued on y.
W_0 = reshape(output.W, nb, ny)'; % W(i,b) expected continuation value on y starting with b.
bgrid = output.bgrid;
ygrid = output.ygrid;
Vd = output.Vd;

Q_10 = zeros(ny, nb); % Store all the values from q_10.
W_10 = zeros(ny, nb); % Store all the values from w_10.


%% Step: Perform one interation given (i,b).

% This step return an updated value for E_m(V^d(y),V^nd(y,b,m)).
% we get policies for d(y,b,m), b(y,b,m), m_grid(y,b,m).
% With this iteration we know that in case the agent lands in (12,1,m)
% Defaults for sure. 

%for i=1:ny
%    for b=1:nb
 
for i=1:ny
    for b=1:nb
            
        q_10 = 0;
        w_10 = 0;
        dp_grid = [];
        bp_grid = [];
        m_grid = [];
        Exit = false;
        
        C = zeros(params.nb,1);
        W = zeros(params.nb,1);
        
        for x=1:nb
            C(x) = ygrid(i) - Q_0(i,x) * (bgrid(x) - (1-lambda) * bgrid(b)) + (lambda + z * (1-lambda)) * bgrid(b);
            W(x) = W_0(i,x);
        end
        
        save("Wlist.mat", "W")
        save("clist.mat", "C")
        % Create X_tilde:
        
        X_tilde = [];
        for x=1:nb
            if C(x) + m_bar > 0
                X_tilde(end+1)=x;
            end
        end
        
        % If X_tilda is empty or not:
        if isempty(X_tilde) % X_tilda is empty.
            Exit = true;
            m_grid = [-params.m_bar, params.m_bar];
            dp_grid = 1;
            bp_grid = nan;
            q_10 = 0;
            w_10 = output.Vd(i);
        end
        
        if (~Exit) % X_tilda(m_bar) is not emtpy:
            % Compute the optimal policy:
            temp1 = -inf;
            x_1 = nan;
            %A = [];
            for x=1:nb
                if C(x) + m_bar > 0 
                    temp2 = V_x(gamma, C(x)+m_bar, W(x));
                    disp(temp2)
                    if temp2 > temp1
                        temp1 = temp2;
                        x_1 = x;
                    elseif temp2 == temp1
                        if C(x) <= C(x_1)
                            temp1 = temp2;
                            x_1 = x;
                        end
                    end
                end
                %A(end+1)= V_x(gamma, C(x)+m_bar, W(x));
            end
            m_i = m_bar;
            x_i = x_1;
        end
        
        while (~Exit)
            % Start by checking if it is better to default:
            if V_x(gamma, C(x_i) + m_i, W(x_i)) <= Vd(i)
                dp_grid = [1, dp_grid];
                m_grid = [-m_bar, m_i, m_grid];
                bp_grid = [nan, bp_grid];
                w_10 = w_10 + cdf_tr(m_i, -m_bar, m_bar, sigma_m) * Vd(i);
                Exit = true;
                break
            else
                X_tilde = []; % {x: c_x > c_i}
                for x=1:nb
                    if C(x) > C(x_i)
                        X_tilde(end+1)=x;
                    end
                end
                if isempty(X_tilde) % There is nothing to replace the policy with.
                   % FUNCTION_LAST_UPDATE:
                   if C(x_i)-m_bar>=10e-6 %1.
                       if V_x(gamma, C(x_i) - m_bar, W(x_i)) >= Vd(i)
                           m_grid = [-m_bar, m_i, m_grid];
                           dp_grid = [0, dp_grid];
                           bp_grid = [x_i, bp_grid];
                           Aux = quadgk(@(m)((V_x(gamma, C(x_i) + m, W(x_i))).*normpdf(m,0,sigma_m)),-m_bar,m_i);
                           w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                           q_10 = q_10 + cdf_tr(m_i, -m_bar, m_bar, sigma_m) * (lambda + (1-lambda) * (Q_0(i,x_i) + z)) * (1/(1+r));
                           Exit = true;
                           break
                       else
                           m_critical = ((Vd(i) - W(x_i)) * (1-gamma))^(1/(1-gamma)) - C(x_i);
                           m_grid = [-m_bar, m_critical, m_i, m_grid];
                           dp_grid = [1, 0, dp_grid];
                           bp_grid = [nan, x_i, bp_grid];
                           Aux = quadgk(@(m)((V_x(gamma, C(x_i)+m, W(x_i))).*normpdf(m,0,sigma_m)), m_critical, m_i);
                           w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                           w_10 = w_10 + Vd(i) * cdf_tr(m_critical, -m_bar, m_bar, sigma_m);
                           q_10 = q_10 + (cdf_tr(m_i, -m_bar, m_bar, sigma_m)-cdf_tr(m_critical, -m_bar, m_bar, sigma_m)) * (lambda + (1-lambda) * (Q_0(i, x_i) + z)) * (1/(1+r));
                           Exit = true;
                           break
                       end
                   else %2.
                       m_critical = ((Vd(i) - W(x_i)) * (1-gamma))^(1/(1-gamma)) - C(x_i);
                       m_grid = [-m_bar, m_critical, m_i, m_grid];
                       dp_grid = [1, 0, dp_grid];
                       bp_grid = [nan, x_i, bp_grid];
                       Aux = quadgk(@(m)((V_x(gamma, C(x_i)+m, W(x_i))).*normpdf(m,0,sigma_m)),m_critical, m_i);
                       w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                       w_10 = w_10 + Vd(i) * cdf_tr(m_critical, -m_bar, m_bar, sigma_m);
                       q_10 = q_10 + (cdf_tr(m_i, -m_bar, m_bar, sigma_m)-cdf_tr(m_critical, -m_bar, m_bar, sigma_m)) * (lambda + (1-lambda) * (Q_0(i,x_i) + z)) * (1/(1+r));
                       Exit = true;
                       break
                   end
                else % There is something to replace the policy with.
                   % FUNCTION_UPDATE:
                   % First, we find m_hat:
                   m_hat = -inf;
                   x_hat = nan;
                   for x=1:nb
                        if C(x) > C(x_i)
                            m_candidate = Bisect(gamma, m_i, C(x_i), W(x_i), C(x), W(x));
                            if m_candidate > m_hat
                                m_hat = m_candidate;
                                x_hat = x;
                            elseif m_candidate == m_hat
                                if C(x)>=C(x_hat)
                                    m_hat = m_candidate;
                                    x_hat = x;
                                end
                            end
                        end
                   end
        
                   if m_hat > -m_bar
                        if V_x(gamma, C(x_i) + m_hat, W(x_i))>= Vd(i)
                               m_grid = [m_i, m_grid];
                               dp_grid = [0, dp_grid];
                               bp_grid = [x_i, bp_grid];
                               Aux = quadgk(@(m)((V_x(gamma, C(x_i)+m, W(x_i))).*normpdf(m,0,sigma_m)), m_hat, m_i);
                               w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                               q_10 = q_10 + (cdf_tr(m_i, -m_bar, m_bar, sigma_m)-cdf_tr(m_hat, -m_bar, m_bar, sigma_m)) * (lambda + (1-lambda) * (Q_0(i, x_i) + z)) * (1/(1+r));
                               x_i = x_hat;
                               m_i = m_hat;
                               Exit = false;
                        else
                               m_critical = ((Vd(i) - W(x_i)) * (1-gamma))^(1/(1-gamma)) - C(x_i);
                               m_grid = [-m_bar, m_critical, m_i, m_grid];
                               dp_grid = [1, 0, dp_grid];
                               bp_grid = [nan, x_i, bp_grid];
                               Aux = quadgk(@(m)((V_x(gamma, C(x_i)+m, W(x_i))).*normpdf(m,0,sigma_m)),m_critical, m_i);
                               w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                               w_10 = w_10 + Vd(i) * cdf_tr(m_critical, -m_bar, m_bar, sigma_m);
                               q_10 = q_10 + (cdf_tr(m_i, -m_bar, m_bar, sigma_m)-cdf_tr(m_critical, -m_bar, m_bar, sigma_m)) * (lambda + (1-lambda) * (Q_0(i, x_i) + z)) * (1/(1+r));
                               Exit = true;
                               break
                         end
                   else
                       if V_x(gamma, C(x_i)-m_bar, W(x_i))>= Vd(i)
                           m_grid = [-m_bar, m_i, m_grid];
                           dp_grid = [0, dp_grid];
                           bp_grid = [x_i, bp_grid];
                           Aux = quadgk(@(m)((V_x(gamma, C(x_i)+ m, W(x_i))).*normpdf(m,0,sigma_m)),-m_bar,m_i);
                           w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                           q_10 = q_10 + cdf_tr(m_i, -m_bar, m_bar, sigma_m) * (lambda + (1-lambda) * (Q_0(i,x_i)+z)) * (1/(1+r));
                           Exit = true;
                           break
                       else
                           m_critical = ((Vd(i) - W(x_i)) * (1-gamma))^(1/(1-gamma)) - C(x_i);
                           m_grid = [-m_bar, m_critical, m_i, m_grid];
                           dp_grid = [1, 0, dp_grid];
                           bp_grid = [nan, x_i, bp_grid];
                           Aux = quadgk(@(m)((V_x(gamma, C(x_i) + m, W(x_i))).*normpdf(m,0,sigma_m)),m_critical, m_i);
                           w_10 = w_10 + Aux / (normcdf(m_bar, 0, sigma_m)-normcdf(-m_bar, 0, sigma_m));
                           w_10 = w_10 + Vd(i) * cdf_tr(m_critical, -m_bar, m_bar, sigma_m);
                           q_10 = q_10 + (cdf_tr(m_i, -m_bar, m_bar, sigma_m)-cdf_tr(m_critical, -m_bar, m_bar, sigma_m)) * (lambda + (1-lambda) * (Q_0(i,x_i)+z)) * (1/(1+r));
                           Exit = true;
                           break
                       end
                   end
                end
            end
        end
        Q_10(i,b)=q_10;
        W_10(i,b)=w_10;
    end
end


%% Take the expectations over y:
Q_1 = zeros(ny,nb);
W_1 = zeros(ny,nb);

for i=1:ny
    for b=1:nb
        for i_prime=1:ny
            Q_1(i,b) = Q_1(i,b)+ P(i,i_prime)*Q_10(i_prime,b);s
            W_1(i,b) = W_1(i,b)+ beta * P(i,i_prime)*W_10(i_prime,b);
        end
    end
end

disp(norm(Q_1 - Q_0, 'inf'));
disp(norm(W_1 - W_0, 'inf'));
Q2 = Q_1 - Q_0;
W_2 = W_1 - W_0;
%% Functions:

% Value of repayment under policy x:
function Value = V_x(gamma, c_x, W_x)
    Value = c_x.^(1-gamma)/(1-gamma) + W_x;
end

% Root finding for c_candidate: 
function m_indf = Bisect(gamma, m_i, c_i, W_i, c_x, W_x)
    % need to handle both cases 3 and 4
    mlb = -c_i + 10e-8;
    mub = m_i;
    err_mub = V_x(gamma, c_i + mub, W_i) - V_x(gamma, c_x + mub, W_x); % >0 by assumption.
    bContinue = true;
    while bContinue
        m_switch = (mlb+mub)/2;
        err_mswitch = V_x(gamma, c_i + m_switch, W_i) - V_x(gamma, c_x + m_switch, W_x);
        if sign(err_mub)*sign(err_mswitch)<=0 
            mlb = m_switch;
        else
            mub = m_switch;
            err_mub = err_mswitch;
        end
        err = mub-mlb;
        if err< 10e-8
            bContinue = false;
        end
    end
    m_indf = (mlb + mub)/2;

end

% Cumulative probabilities over [-m_bar, x] of the truncated normal:
function cum_prob = cdf_tr(x, a, b, sigma)
    denominator = (normcdf(b, 0, sigma) - normcdf(a, 0, sigma));
    if x>=a && x<=b 
        cum_prob = (normcdf(x, 0, sigma)-normcdf(a, 0, sigma))/denominator;
    else
        cum_prob = nan;
    end
end







