
% set number of force curve and spans of model parameters
n_force_curve = 1; % 1..6
E1_span = [7e2; 4e-0; 2.8e3; 4e-0; 4e-0; 4e-0];
tau1_span = [5e-1; 1e0; 1e-2; 1e0; 1e0; 1e0];
theta_span = [1e0; 1e-3; 1e-2; 1e-3; 1e-3; 1e-3];
eps_loading_span = [1e2; 1; .005; 1; 1; 1];

% model parameters
lambda = 1e0;
gamma = 1e0;
E0_loading = 3e2; % sigma_dwell_frac(end) / eps_dwell(end)
% E0_dwell = E0_loading; % 
E1 = E1_span(n_force_curve);
tau_1 = tau1_span(n_force_curve);
nu1 = tau_1 * E1;
R = 5e3;

% time
time_step = .0013;
time = (time_step:time_step:10)';
time_loading = time(time <= 1);
time_dwell = time(time >= 1);

% deformation|identation
eps_loading = eps_loading_span(n_force_curve);
eps_dwell = eps_loading_span(n_force_curve);
deps_dt = eps_loading * time_step;
deps_dt_dwell = 0;

% for initial conditions
F_loading = 4/3 * sqrt(R) * (E1 * exp(- time_loading / tau_1) + E0_loading) .* (eps_loading * time_loading).^(3/2);
sigma_loading_rel_fun = F_loading ./ (eps_loading * time_loading).^(3/2) .* (eps_loading * time_loading);

% effective temperature factor
theta = theta_span(n_force_curve);

% approximaxion model
[~,fit_model] = get_xi_approximation;

% loading
initial_conditions = sigma_loading_rel_fun(1);
[~,sigma_loading] = ode15s(@(t, sigma)...
    get_constitutive_equation_rhs(t, sigma, eps_loading, deps_dt, E0_loading, E1, nu1, theta, lambda, gamma, fit_model), time_loading, initial_conditions);
sigma_o_loading = sigma_loading - E0_loading * (eps_loading * time_loading);
eps_o_loading = eps_loading * time_loading - sigma_o_loading / E1;

% dwell
% initial_conditions = sigma_dwell_frac(1);
% [~,sigma_dwell] = ode15s(@(t, sigma)...
%     get_constitutive_equation_rhs(t, sigma, eps_dwell, deps_dt_dwell, E0_dwell, E1, nu1, theta, lambda, gamma, fit_model), time_dwell, initial_conditions);
% sigma_o_dwell = sigma_dwell - E0_dwell * eps_dwell;
% eps_o_dwell = eps_dwell - sigma_o_dwell / E1;

%%
% visualization
figure(1);hold on;
plot(time_loading,sigma_loading_rel_fun) % ,'k--'
% plot(time_loading,sigma_loading,'k')
% plot(time_loading,sigma_o_loading,'r')
% plot(time_loading,E0_loading * (eps_loading * time_loading),'b')
xlabel('t')
ylabel('\sigma');

% figure(2);hold on;
% plot(time_dwell,sigma_dwell_frac,'k--')
% plot(time_dwell,sigma_dwell,'k')
% % plot(time_dwell,E0_loading * eps_dwell * time_dwell,'r')
% % plot(time_dwell,sigma_dwell - E0_loading * eps_dwell,'r')
% xlabel('t')
% ylabel('\sigma');

% figure(3);hold on;
% plot(sigma_loading - sigma_o_loading, eps_loading * time_loading,'b');
% plot(sigma_o_loading, eps_o_loading,'r');
% xlabel('\sigma_o')
% ylabel([char(949) '_o']);

% figure(4);hold on;
% plot(sigma_o_dwell);

% figure(5);hold on;
% % plot(sigma_dwell_frac - sigma_dwell_zener);
% % plot(sigma_dwell_frac - sigma_dwell);
% plot(sigma_dwell - sigma_dwell_zener)
  
%% Вспомогательные функции
function [dsigma_dt] = get_constitutive_equation_rhs(t, sigma, deps, deps_dt, E0, E1, nu1, theta, lambda, gamma, fit_model)
% Функция, определяющая систему дифф. уравнений для моделирования чистого
% сдвига.

if nargin < 11
  [~,fit_model] = get_xi_approximation;
end

if deps_dt == 0
  eps = deps / t;
else
  eps = deps;
end

eps_o = eps * t - (sigma - E0 * eps * t) / E1;
dF_deps_o = ...
  2/3 * theta / gamma * feval(fit_model,eps_o) - lambda / gamma * eps_o;
  
dsigma_dt = zeros(size(sigma));
dsigma_dt(1) = ...
  (E0 + E1) * deps_dt ...
  + (E0*E1) / nu1 * (eps * t) ...
  - E1 / nu1 * sigma ...
  ; % + E1 / nu1 * dF_deps_o

end