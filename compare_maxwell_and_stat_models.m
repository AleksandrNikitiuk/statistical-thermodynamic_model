
time = logspace(-3, 3, 2^12);
initial_sigma = 3e0; % 6.0001e4
relaxation_times = logspace(-2,2,5);

% normalize_sigma = zeros(numel(time),numel(relaxation_times));
% for relaxation_time = relaxation_times
% 
%   [~,sigma] = ode15s(@(t, sigma)...
%       sode(t, sigma, relaxation_time), time, initial_sigma);
% 
%   normalize_sigma(:,relaxation_time == relaxation_times) = (sigma) / max(sigma);
% 
%   figure(1);hold on;
%   plot(time(normalize_sigma(:,relaxation_time == relaxation_times) > 1e-5),normalize_sigma(normalize_sigma(:,relaxation_time == relaxation_times) > 1e-5,relaxation_time == relaxation_times),'k')
% end
% 
% figure(1);hold on;
% plot(time,sum(normalize_sigma,2) / max(sum(normalize_sigma,2)),'b')
% % plot(time,mean(normalize_sigma,2) / max(mean(normalize_sigma,2)),'b')
% plot(time,time.^-.16 / max(time.^-.16),'r')
% xlabel('t')
% ylabel('\sigma');

E0 = 1e3; % 5.99999e7
E1 = 1e3; % 1e3
tau_1 = 1e-3;
k_boltzmann = 1.380649e-23;
temperature = 3.1e2;
theta = k_boltzmann * temperature / 2;
nu = tau_1 * E1;
lambda_stat = .4e1; % .86e1 .7 .6
gamma_stat = 2.468083232289583e-21;

[~,fit_model] = get_xi_approximation;

strain = 2.6e-3;
pd = makedist('Exponential',1.5e-2);
sigma_model = zeros(numel(time),2^10);
for i = 1:2^10
  chi = 2e-2 + random(pd);
  [~,sigma_model(:,i)] = ode15s(@(t, sigma)...
      get_rhs_expr(t,sigma,strain,E0,E1,nu,theta,lambda_stat,gamma_stat,fit_model,chi), time, initial_sigma);
% disp(i)
end
% normalize_sigma_model = sigma_model / max(sigma_model);

figure(2);hold on
plot((time),(mean(sigma_model,2)) * 1e1)
% plot(log10(time(log10(time) > -1 & log10(time) < 1)),log10(mean(sigma_model(log10(time) > -1 & log10(time) < 1,:),2)) * 1e1)
xlabel('t') %  ['{\it' char(949) '_o}']
ylabel('\sigma');

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


  
  %% Вспомогательные функции
function [dsigma_dt] = sode(t, sigma, relaxation_time)
% Функция, определяющая систему дифф. уравнений для моделирования чистого
% сдвига.

dsigma_dt = zeros(size(sigma));

dsigma_dt(1) = - sigma / relaxation_time;

end

function [rhs_expr] = get_rhs_expr(t,varibles,indentation_rate,E0,E1,nu,theta,lambda,gamma,fit_model,chi)
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

strain = indentation_rate;

rhs_expr = zeros(size(varibles));
eps_o = strain - ( varibles(1) - E0 * strain ) / E1;
ksi = feval(fit_model,eps_o);

orientation_deformation = unique([-.3:.01:.99 eps_o]);
numbers_orientation_deformation = 1:numel(orientation_deformation);
effective_field_integral_orientation_case = ...
  get_effective_field_integral_orientation_case(orientation_deformation);
F_orientation = ...
  effective_field_integral_orientation_case ...
  - (varibles(1) - E0 * strain) .* orientation_deformation ... % gamma / theta * 
  - (chi)^-1 * orientation_deformation.^2 / 2;
F_orientation = F_orientation - mean(F_orientation);
% if eps_o > .8 || eps_o <-.3
  % F_orientation_barrier = F_orientation_barrier;% + get_energy_barrier(F_orientation,numbers_orientation_deformation(orientation_deformation == eps_o));
% else
  F_orientation_barrier = get_energy_barrier(F_orientation,numbers_orientation_deformation(orientation_deformation == eps_o));
% end

dF_deps_o = 2/3 * ksi - (chi)^-1 * eps_o; % theta / gamma * 
rhs_expr(1) = ...
  E1 / (nu * exp(F_orientation_barrier) ) ... %  
  * ( E0 * strain ...
  - varibles(1) + dF_deps_o ...
  );


% figure(1);hold on;
% plot(orientation_deformation,F_orientation + abs(min(F_orientation)))
% plot(eps_o,F_orientation(orientation_deformation == eps_o) + abs(min(F_orientation)),'o')

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


