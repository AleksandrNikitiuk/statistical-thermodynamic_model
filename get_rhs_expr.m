function [rhs_expr] = get_rhs_expr(t,varibles,eps,deps_dt,E1,E2,Lo,nu,theta,lambda,gamma)
% Функция для расчета правой части дифференциального уравнения
%     Пример использования:
%     time_step = .013;
%     time = (time_step:time_step:10)';
%     time_loading = time(time <= 1);
%     time_dwell = time(time > 1);
% 
%     sigma_0 = 0.;
%     initial_conditions = sigma_0;
% 
%     E1 = 300;
%     E2 = 100;
%     Lo = 1;
%     nu = 100;
% 
%     eps = .005;
%     deps_dt = eps * time_step;
%     deps_dt_dwell = 0;
% 
%     theta = [.2];
%     line_style = {'-','--',':'};
% 
%     for i = 1:length(theta)
%       % loading
%       [~,Sigma] = ode15s(@(t,sigma) get_rhs_expr(t,sigma,eps,deps_dt,E1,E2,Lo,nu,theta(i)),time_loading,initial_conditions);
% 
%       % dwell
%       initial_conditions = Sigma(end);
%       [~,Sigma_dwell] = ode15s(@(t,sigma) get_rhs_expr(t,sigma,eps,deps_dt_dwell,E1,E2,Lo,nu,theta(i)),time_dwell,initial_conditions);
%       Sigma = [Sigma; Sigma_dwell];
% 
%       eps_time = [eps*time_loading; eps * ones(length(time_dwell),1)];
% 
%     %   figure(1);hold on
%     %   plot(eps_time,Sigma,'k','LineWidth',2,'LineStyle',line_style{i})
%     %   xlabel(['{\it' char(949) '}']);
%     %   ylabel('\sigma');
%     %   set_figure;
%     % 
%     %   figure(2);hold on
%     %   plot(time,Sigma ./ eps_time,'k','LineWidth',2,'LineStyle',line_style{i})
%     %   xlabel('time');
%     %   ylabel('E');
%     %   set_figure;
% 
%       figure(3);hold on
%       plot(time,(Sigma ./ eps_time) .* eps_time.^(3/2),'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('time');
%       ylabel('F');
%       set_figure;
%     end
%     % for i = 1:3
%     %   figure(i);legend(compose('{\\it\\Theta}=%.2f',theta));
%     % end


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

[~,fit_model] = get_xi_approximation;
dF_deps = ...
  -4/9 * theta / lambda * feval(fit_model,eps) + gamma / lambda * eps;

ab = fit_model.a * fit_model.b;
cd = fit_model.c * fit_model.d;

rhs_expr(1) = 3/2 * (...
  (E1 + E2 + E2 * Lo / nu ...
    + fit_model.e ...                       %
    + ab * exp(fit_model.b * eps * t) ...   % d/dt(dF/deps)
    + cd * exp(fit_model.d * eps * t) ...   %
  ) * deps_dt ...
  + (E1 + E2) / nu * eps * t ...
  - E2 / nu * ( 2/3 * varibles(1) + dF_deps ) ...
  );
end

