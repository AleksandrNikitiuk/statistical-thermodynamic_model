
time = logspace(-3, 5, 2^12);
initial_strain = -1e-4; % .0e0 .2e0
relaxation_times = logspace(-2,2,5);

E0 = 1e3; % 1e3 1e-1 1e1
E1 = 1e6; % 1e3 1e-1 1e3
tau_1 = 2e-3; % 1e2 1e5
k_boltzmann = 1.380649e-23;
temperature = 3.1e2;
theta = 2.3e-16; % 1 2.3 k_boltzmann * temperature / 2
nu = tau_1 * E1;
lambda_stat = E1; % .4e1 .86e1 .7 .6
gamma_stat = 2e-21; % 2e-21 5e-21
N_stat = gamma_stat^-1;

[~,fit_model] = get_xi_approximation;

stress = -1e-9 * gamma_stat / theta; % 1e3 1e0 1e-2 1e-3 1e-9

chi = theta / ( lambda_stat * gamma_stat);
E1 = E1 * gamma_stat / theta;
E0 = E0 * gamma_stat / theta;
tau_1 = tau_1 * (sqrt(gamma_stat^(1/3) * lambda_stat / 1e-16 ))^1;

pd = makedist('Exponential',chi*1e-1);
n_realization = 2^0;
strain_model = zeros(numel(time),n_realization);
for i = 1:n_realization
  chi_random = chi; %  + random(pd)
  [~,strain_model(:,i)] = ode15s(@(t, strain)...
      get_rhs_expr(t,strain,stress,E0,E1,nu,theta,lambda_stat,gamma_stat,fit_model,chi_random), time, initial_strain);
disp(i)
end
% normalize_sigma_model = sigma_model / max(sigma_model);

% x = linspace(0,chi*1e-1,2^10);
% figure(3);
% plot(pdf(pd,chi + x),chi + x)

figure(2);hold on
plot((time),.01*((mean(strain_model,2))))
% plot(log10(time(log10(time) > -1 & log10(time) < 1)),log10(mean(sigma_model(log10(time) > -1 & log10(time) < 1,:),2)) * 1e1)
xlabel('t') %  ['{\it' char(949) '_o}']
ylabel(['{\it' char(949) '}']);

% chi_orientation = .4; % .06 .1 .2
% orientation_deformation = (0:.01:.4)';
% [~,fit_model] = get_xi_approximation;
% initial_sigma = 7;
% sigma_orientation = initial_sigma - 1 * (feval(fit_model,orientation_deformation) - chi_orientation^-1 * orientation_deformation);
% 
% figure(1);hold on;
% plot((logspace(-3,3,numel(sigma_orientation))),sigma_orientation,'LineWidth',3);
% ylabel('{\it\sigma}');
% xlabel(['{\itt}']);

% dimensional parameters
% E0 * (chi * lambda_stat * gamma_stat)  / gamma_stat
% E1 * (chi * lambda_stat * gamma_stat)  / gamma_stat
% tau_1 * (sqrt(gamma_stat^(1/3) * lambda_stat / 1e-22 ))^-1


  
  %% Вспомогательные функции
function [rhs_expr] = get_rhs_expr(t,varibles,stress,E0,E1,nu,theta,lambda,gamma,fit_model,chi)
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

% stress = stress_rate;

eps_o = (varibles(1) - ( stress - E0 * varibles(1) ) / E1);
ksi = feval(fit_model,eps_o) - 0.295369501433937;

orientation_deformation = unique([-.29:.01:.79 eps_o]);
numbers_orientation_deformation = 1:numel(orientation_deformation);
effective_field_integral_orientation_case = ...
  get_effective_field_integral_orientation_case(orientation_deformation);
F_orientation = ...
  effective_field_integral_orientation_case ...
  - (stress - E0 * varibles(1)) .* orientation_deformation ... % gamma / theta * 
  - (chi)^-1 * orientation_deformation.^2 / 2; %  - orientation_deformation.^-1
F_orientation = F_orientation - mean(F_orientation);
% if eps_o > .8 || eps_o <-.3
  % F_orientation_barrier = F_orientation_barrier;% + get_energy_barrier(F_orientation,numbers_orientation_deformation(orientation_deformation == eps_o));
% else
  F_orientation_barrier = get_energy_barrier(F_orientation,numbers_orientation_deformation(orientation_deformation == eps_o));
% end

dF_deps_o = 2/3 * ksi - (chi)^-1 * eps_o; % theta / gamma * 
rhs_expr = zeros(size(varibles));
rhs_expr(1) = E1 / (nu * ( E0 + E1 ))... %  * exp(F_orientation_barrier)
  * (...
  stress - E0 * varibles(1) - dF_deps_o ...
  );


figure(1);hold on;
plot(orientation_deformation,F_orientation + abs(min(F_orientation)),'k','LineWidth',3)
plot(eps_o,F_orientation(orientation_deformation == eps_o) + abs(min(F_orientation)),'ob','MarkerSize',24,'LineWidth',3)
% xlim([-.15e-2 .15e-2])
% ylim([0 .08])
% plot((t),abs(varibles(1)),'o')

% sigma_reaction = (feval(fit_model,orientation_deformation)' - (theta / (lambda * gamma))^-1 * orientation_deformation);
% plot(log10(sigma_reaction + abs(min(sigma_reaction))),log10(orientation_deformation + abs(min(orientation_deformation))))
% plot(log10(sigma_reaction(orientation_deformation == eps_o) + abs(min(sigma_reaction))),log10(eps_o + abs(min(orientation_deformation))),'o')
% xlim(log10([-3 2] + abs(min(sigma_reaction))))
end

function [effective_field_derivative] = get_effective_field_integral_orientation_case(deformation)

effective_field_derivative = ...
  - 1.7036 * log(0.993712 - deformation)...
 - 114161 * log(8822.34 - deformation) ...
 - 0.911178 * log(0.498379 + deformation) ...
 - 36.3414 * log(2.83803 + deformation) ...
  ;
effective_field_derivative = real(effective_field_derivative);

end

function [energy_barrier] = get_energy_barrier(free_energy,ind)

% free_energy = free_energy(ind:end) + abs(min(free_energy(ind:end)));
% inds = 1:numel(free_energy);
% 
% local_minima_inds = inds(islocalmin(free_energy));
% local_maxima_inds = inds(islocalmax(free_energy));
% 
% if isempty(local_minima_inds) && isempty(local_maxima_inds)
%   % energy_barrier = abs(free_energy(1) - free_energy(end));
%   if free_energy(1) < free_energy(end)
%     energy_barrier = free_energy(end) - free_energy(1);
%   else
%     energy_barrier = 0;
%   end
%   return;
% end
% 
% if isempty(local_minima_inds) && ~isempty(local_maxima_inds)
%   [~,local_minima_inds] = min(free_energy);
% else
%   [~,local_maxima_inds] = max(free_energy);
% end
% 
% if local_minima_inds(1) < local_maxima_inds(1)
%   energy_barrier = free_energy(local_minima_inds(1)) + free_energy(1);
% else
%   energy_barrier = free_energy(local_maxima_inds(1)) - free_energy(1);
% end

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


