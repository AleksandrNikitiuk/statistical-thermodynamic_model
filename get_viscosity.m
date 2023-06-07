function [viscosity] = get_viscosity(initial_viscosity,eps_o,sigma_o,theta,lambda,gamma)
% Функция для расчета коэффициента вязкости.
%     Пример использования:
%     % time
%     time_step = .013;
%     time = (time_step:time_step:10)';
%     time_loading = time(time <= 1);
%     time_dwell = time(time > 1);
% 
%     % initial conditions
%     eps_o_0 = 0.001;
% 
%     % model parameters
%     tau_1 = 1;
%     E1 = 1000;
%     Lo = 1 / (tau_1 * E1);
% 
%     lambda = 1;
%     gamma = 1/E1;
% 
%     % deformation
%     eps = .5;
%     deps_dt = eps * time_step;
% 
%     % temperature
%     theta = [.2 .3 .4];
%     line_style = {'-','--',':'};
% 
%     for i = 1:length(theta)
% 
%       initial_conditions = eps_o_0;
%       [~,eps_o_loading] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps,deps_dt,E1,0,Lo,0,theta(i),lambda,gamma),time,initial_conditions);
%       sigma_o_loading =  E1 * (eps * time - eps_o_loading);
% 
%       figure(1);hold on
%       yyaxis left
%       plot(time,eps_o_loading,'b','LineWidth',2,'LineStyle',line_style{i})
%       ylabel([char(949) '_o']);
%       yyaxis right
%       plot(time,sigma_o_loading,'r','LineWidth',2,'LineStyle',line_style{i})
%       ylabel('\sigma_o');
%       xlabel('time');
%       set_figure;
% 
%       figure(2);hold on
%       plot(eps_o_loading,sigma_o_loading,'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel([char(949) '_o']);
%       ylabel('\sigma_o');
%       set_figure;
%     end


if nargin < 2
  viscosity = initial_viscosity;
  return
end

if nargin < 4
  theta = .2;
end

if nargin < 6
  lambda = 1;
  gamma = 1;
end

viscosity = initial_viscosity ...
  * exp( (lambda * eps_o.^2 - gamma * sigma_o) / theta); % 

end

