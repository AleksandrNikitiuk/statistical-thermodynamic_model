function [rhs_expr] = get_rhs_expr(t,varibles,eps,deps_dt,E1,E2,Lo,nu,theta,lambda,gamma)
% Функция для расчета правой части дифференциального уравнения
%     Пример использования:
%     % time
%     time_step = .013;
%     time = (time_step:time_step:10)';
%     time_loading = time(time <= 1);
%     time_dwell = time(time > 1);
% 
%     % initial conditions 
%     sigma_0 = 0;
%     eps_o_0 = 0;
% 
%     % model parameters
%     E0 = 2000;
% 
%     tau_1 = .5;
%     E1 = 1000;
%     Lo = 1 / (tau_1 * E1);
% 
%     tau_2 = .5;
%     E2 = 1000;
%     nu = tau_2 * E2;
% 
%     lambda = 1;
%     gamma = 1;
% 
%     % deformation
%     eps = .005;
%     deps_dt = eps * time_step;
%     deps_dt_dwell = 0;
% 
%     % temperature
%     theta = [.2 .3 .4];
%     line_style = {'-','--',':'};
% 
%     for i = 1:length(theta)
%       % loading regime
%       sigma_e_loading = E0 * eps * time_loading;
% 
%       initial_conditions = eps_o_0;
%       [~,eps_o_loading] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps,deps_dt,E1,0,Lo,0,theta(i),lambda,gamma),time_loading,initial_conditions);
%       sigma_o_loading =  E1 * (eps * time_loading - eps_o_loading);
% 
%       initial_conditions = sigma_0;
%       [~,sigma_m_loading] = ode15s(@(t,sigma_m) get_rhs_expr(t,sigma_m,eps,deps_dt,0,E2,0,nu,theta(i),lambda,gamma),time_loading,initial_conditions);
% 
%       % dwell regime
%       sigma_e_dwell = E0 * eps * time_loading(end) * ones(size(time_dwell));
% 
%       initial_conditions = eps_o_loading(end);
%       [~,eps_o_dwell] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps,deps_dt_dwell,E1,0,Lo,0,theta(i),lambda,gamma),time_dwell,initial_conditions);
%       sigma_o_dwell =  E1 * (eps * time_loading(end) - eps_o_dwell);
% 
%       initial_conditions = sigma_m_loading(end);
%       [~,sigma_m_dwell] = ode15s(@(t,sigma_m) get_rhs_expr(t,sigma_m,eps,deps_dt_dwell,0,E2,0,nu,theta(i),lambda,gamma),time_dwell,initial_conditions);
% 
%       % stress
%       sigma_e = [sigma_e_loading; sigma_e_dwell];
%       eps_o = [eps_o_loading; eps_o_dwell];
%       sigma_o = [sigma_o_loading; sigma_o_dwell];
%       sigma_m = [sigma_m_loading; sigma_m_dwell];
%       sigma = sigma_e + sigma_o + sigma_m;
% 
%       % deformation
%       eps_time = [eps*time_loading; eps * ones(size(time_dwell))];
% 
%       % Young modulus
%       E = sigma ./ eps_time;
%       E_e = sigma_e ./ eps_time;
%       E_o = sigma_o ./ eps_time;
%       E_m = sigma_m ./ eps_time;
% 
%       % Force
%       F = E .* eps_time.^(3/2);
%       F_e = E_e .* eps_time.^(3/2);
%       F_o = E_o .* eps_time.^(3/2);
%       F_m = E_m .* eps_time.^(3/2);
% 
%       % Visualization
%       figure(1);hold on
%       plot(eps_time,sigma,'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel(['{\it' char(949) '}']);
%       ylabel('\sigma');
%       set_figure;
% 
%       figure(2);hold on
%       plot(time,E,'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('time');
%       ylabel('E');
%       set_figure;
% 
%       figure(3);hold on
%       plot(time,F,'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('time');
%       ylabel('F');
%       set_figure;
%     end
% 
%     % Add legend
%     for i = 1:3
%       figure(i);legend(compose('{\\it\\Theta}=%.2f',theta));
%     end



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
  sigma = E1 * (eps * t - varibles(1));
  
  [~,fit_model] = get_xi_approximation;
  dF_deps_o = ...
    -2/3 * theta / gamma * feval(fit_model,varibles(1)) + lambda / gamma * varibles(1);

  rhs_expr(1) = 1/get_viscosity(1/Lo) * ( sigma + dF_deps_o ); % ,varibles(1)
else
  rhs_expr(1) = E2 * deps_dt - E2/nu * varibles(1);
end

end

