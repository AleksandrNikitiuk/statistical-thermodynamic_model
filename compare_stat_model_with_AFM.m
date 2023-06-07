%% Сравнение стат. модели и результатов АСМ.


% set number of force curve and spans of model parameters
n_force_curve = 1; % 1..6
E1_span = [3.6e-1; 90; 450; 700; 900; 1700]; % [E1] = 10 kPa
tau1_span = [1e-0; 1; 1; 2; 3; 2];
theta_span = [6e-2; .16; .16; .16; .16; .16]; % [theta] = eV
eps_loading_span = [.2; .22; .25; .3; .33; .39];

% load AFM data
n_point = 1;
load("force_ind_time_cell08speed_.mat");
time = time(indentation(:,n_point) >= 0);
time = time - min(time); % s
force = force(indentation(:,n_point) >= 0,n_point); % [force] = nN
indentation = indentation(indentation(:,n_point) >= 0,n_point); % [indentation] = nm

% time
[~,index_max_indentation_time] = max(indentation);
time_loading = time(time <= time(index_max_indentation_time));
time_unloading = time(time > time(index_max_indentation_time));
time_step = time(2) - time(1);

% deformation|identation
eps_loading = indentation(1:index_max_indentation_time) / max(indentation); % [eps] = 1L
eps_unloading = indentation((index_max_indentation_time+1):end) / max(indentation);
deps_dt = eps_loading * time_step; % [deps/dt] = s^-1
deps_dt_unloading = eps_unloading * time_step;

% model parameters
E1 = E1_span(n_force_curve);
tau_1 = tau1_span(n_force_curve);
Lo = 1 / (tau_1 * E1);
lambda = 1e-0;
gamma = 1e-0; % TODO: не должно зависеть от E1, размерность L^-3, но порядки д.б. обратно порпорциональны

% initial conditions
F = force; %  / max(force)
F_loading = F(1:index_max_indentation_time);
F_unloading = F((index_max_indentation_time+1):end);
sigma_0_loading = F_loading(1) / eps_loading(1)^(3/2) * eps_loading(1);
sigma_0_unloading = F_unloading(1) / eps_unloading(1)^(3/2) * eps_unloading(1);
eps_o_0_loading = 0; % eps_loading(1) - sigma_0_loading / E1
eps_o_0_unloading = eps_unloading(1) - sigma_0_unloading / E1;

% model parameters
E0_loading = 1e-4; % [E0] = 10 kPa
E0_unloading = E0_loading; % F_unloading(1) / eps_unloading(1) * 10^-1
 
% effective temperature factor
theta = theta_span(n_force_curve);

% loading regime
initial_conditions = eps_o_0_loading;
[~,eps_o_loading] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps_loading,deps_dt,time_loading,E1,0,Lo,0,theta,lambda,gamma),time_loading,initial_conditions);
sigma_o_loading =  E1 * (eps_loading - eps_o_loading);

% unloading regime
initial_conditions = eps_o_loading(end);
[~,eps_o_unloading] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps_unloading,deps_dt_unloading,time_unloading,E1,0,Lo,0,theta,lambda,gamma),time_unloading,initial_conditions);
sigma_o_unloading =  E1 * (eps_unloading - eps_o_unloading);

%%
% visualization results
n = 2.1e-1; % 7.8e4
F_loading_stat = ((E0_loading * eps_loading + sigma_o_loading) ./ (eps_loading) .* (eps_loading).^(3/2)) / n;
F_unloading_stat = ((E0_unloading * eps_unloading + sigma_o_unloading) ./ (eps_unloading) .* (eps_unloading).^(3/2)) / n;
F_stat = [F_loading_stat; F_unloading_stat];

figure(1);hold on
% plot(time_loading,sigma_o_loading,'r','LineWidth',2)
plot(...
  (eps_loading),...
  F_loading_stat,... % 
  'k','LineWidth',2)
plot((eps_loading),(F_loading),'k--','LineWidth',2)
xlabel('indentation, nm');
ylabel('{\itF}, nN');
legend({'stat. model','AFM data'},'Location','northwest');
set_figure;

figure(2);hold on;
plot(E0_loading * eps_loading,eps_loading - eps_o_loading,'b')
plot(sigma_o_loading, eps_o_loading,'r')
plot(E0_loading * eps_loading + sigma_o_loading,eps_loading,'k')
xlabel('\sigma');
ylabel('eps');

figure(3);hold on
% plot(time_unloading,sigma_o_unloading,'r','LineWidth',2)
plot(...
  (eps_unloading),...
  F_unloading_stat,... % 
  'k','LineWidth',2)
plot((eps_unloading),(F_unloading),'k--','LineWidth',2)
xlabel('time, s');
ylabel('{\itF}, nN');
legend({'stat. model','AFM data'},'Location','northeast');
set_figure;

figure(4);hold on;
plot(E0_unloading * eps_unloading,eps_unloading - eps_o_unloading,'b')
plot(sigma_o_unloading, eps_o_unloading,'r')
plot(E0_unloading * eps_unloading + sigma_o_unloading,eps_unloading,'k')
xlabel('\sigma');
ylabel('eps');

% figure(5);hold on
% % plot(time_loading,sigma_o_loading,'r','LineWidth',2)
% plot(...
%   (time),...
%   F_stat,... % 
%   'k','LineWidth',2)
% plot((time),(F),'k--','LineWidth',2)
% xlabel('time, s');
% ylabel('{\itF}, nN');
% legend({'stat. model','AFM data'},'Location','northeast');
% set_figure;

figure(6);hold on;
plot(indentation,F_stat,'k','LineWidth',2);
plot((indentation),(F),'k--','LineWidth',2)
xlabel('indentation, nm');
ylabel('{\itF}, nN');
% legend({'stat. model','AFM data'},'Location','northeast');
set_figure;

function [rhs_expr] = get_rhs_expr(t,varibles,eps,deps_dt,time,E1,E2,Lo,nu,theta,lambda,gamma)
% Функция для расчета правой части дифференциального уравнения
if nargin < 11
  lambda = 1.;
  gamma = 1.;
end

if nargin < 9
  theta = .2;
end

if nargin < 8
  nu = 1.;
end

if nargin < 7
  Lo = 1.;
end

if nargin < 6
  E2 = 1.;
end

if nargin < 5
  E1 = 1.;
end

if nargin < 4
  deps_dt = 0.;
end

if deps_dt == 0
  eps = eps / t;
end

rhs_expr = zeros(size(varibles));

if E2 == 0
  sigma = E1 * (get_eps(t,time,eps) - varibles(1));
  
  [~,fit_model] = get_xi_approximation;
  dF_deps_o = ...
    -2/3 * theta / gamma * feval(fit_model,varibles(1)) + lambda / gamma * varibles(1);

  rhs_expr(1) = 1/get_viscosity(1/Lo,varibles(1),sigma,theta,lambda,gamma) * ( sigma + dF_deps_o ); 
else
  rhs_expr(1) = E2 * deps_dt - E2/nu * varibles(1);
end

end

function eps = get_eps(t,time,epsilons)

eps = epsilons(time == t);

if isempty(eps) || length(eps) ~= 1
  fit_model = fit(time',epsilons,'poly8');
  eps = feval(fit_model,t);
  if eps < epsilons(1)
    eps = epsilons(1);
  end
end

end


