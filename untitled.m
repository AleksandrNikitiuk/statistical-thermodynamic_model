E0 = 1e0;
E1 = 1e2; % 1000 1500 1500 1700
tau_1 = 1.2886e-1; % 1e0 1e0


k_boltzmann = 1.380649e-23;                                                 % Дж/К
temperature = 300;                                                          % К
theta = k_boltzmann * temperature / 2;
nu = tau_1 * E1;
lambda_stat = 2e1;
gamma_stat = 1e-20;

% stat. model
[~,fit_model] = get_xi_approximation;

eps_0 = 1e-5;
initial_conditions = eps_0;
loading_rate = 1;

time_for_stat_model = linspace(0,1,2^12)';
time_step = time_for_stat_model(2) - time_for_stat_model(1);
[~,eps] = ode15s(@(t,eps) get_rhs_expr(t,eps,loading_rate,E0,E1,nu,theta,lambda_stat,gamma_stat,fit_model),time_for_stat_model,initial_conditions);

figure(1);hold on;
plot(time_for_stat_model,eps);



function [rhs_expr] = get_rhs_expr(t,varibles,loading_rate,E0,E1,nu,theta,lambda,gamma,fit_model)
% Функция для расчета правой части дифференциального уравнения

if nargin < 10
  [~,fit_model] = get_xi_approximation;
end

if nargin < 9
  gamma = 1.;
end

if nargin < 8
  lambda = 1.;
end

if nargin < 7
  theta = .2;
end

if nargin < 6
  nu = 1.;
end

if nargin < 5
  E1 = 1.;
end

if nargin < 4
  E0 = 1.;
end

if nargin < 3
  loading_rate = 1;
end

rhs_expr = zeros(size(varibles));

% if E2 == 0
  % sigma_t = loading_rate * t;
  sigma_t = (.2 / (pi * 1 * varibles(1).^2)) * loading_rate * t;
  eps_o = varibles(1) - ( sigma_t - E0 * varibles(1) ) / E1;
  ksi = feval(fit_model,eps_o);
  
  orientation_deformation = -.3:.01:.8;                                       % 1
  effective_field_integral_orientation_case = ...
    get_effective_field_integral_orientation_case(orientation_deformation);
  F_orientation = ...
    effective_field_integral_orientation_case ...
    - gamma / theta * (sigma_t - E0 * varibles(1)) .* orientation_deformation ...
    - (theta / (lambda * gamma))^-1 * orientation_deformation.^2 / 2;
  F_orientation = F_orientation - mean(F_orientation);
  F_orientation_barrier = get_energy_barrier(F_orientation);
  
  dF_deps_o = 2/3 * theta / gamma * ksi - lambda * eps_o;

  % rhs_expr(1) = E1 / (nu * 1 * (E0 + E1)) * ...
  %   ( nu * 1 / E1 * loading_rate + sigma_t ...
  %   - E0 * varibles(1) - 0 );
  rhs_expr(1) = ...
    ( nu / E1 * (.2 / (pi * 1 * varibles(1).^2)) * loading_rate ...
    + sigma_t ...
    - E0 * varibles(1) - 0 ) ...
    / ...
    ( (nu * 1 * (E0 + E1)) / E1  + nu / E1 * .4 / (pi * 1 * varibles(1).^3) * loading_rate * t );
% else
%   rhs_expr(1) = E2 * deps_dt - E2/nu * varibles(1);
% end

end

function [effective_field_derivative] = get_effective_field_integral_orientation_case(deformation)

effective_field_derivative = ...
  - 1.7036 * log(0.993712 - deformation)...
 - 114161 * log(8822.34 - deformation) ...
 - 0.911178 * log(0.498379 + deformation) ...
 - 36.3414 * log(2.83803 + deformation) ......
  ;
effective_field_derivative = real(effective_field_derivative);

end

function [energy_barrier] = get_energy_barrier(free_energy)

local_minima = free_energy(islocalmin(free_energy));

if isempty(local_minima)
  local_maxima = free_energy(islocalmax(free_energy));
  if isempty(local_maxima)
    [~,n_maximum] = max(free_energy);
    if n_maximum == 1
      energy_barrier = 0;
    else
      energy_barrier = free_energy(n_maximum) - min(free_energy);
    end
  else
    energy_barrier = local_maxima(1) - min(free_energy(free_energy < local_maxima(1)));
  end
else
  local_maxima = free_energy(islocalmax(free_energy));
  if isempty(local_maxima)
    free_energy = free_energy(find(free_energy == local_minima(1)):end);
    [~,n_maximum] = max(free_energy);
    if free_energy(n_maximum) == local_minima(1)
      energy_barrier = 0;
    else
      energy_barrier = free_energy(n_maximum) - local_minima(1);
    end
  else
    energy_barrier = local_maxima(1) - local_minima(1);
  end
end

end