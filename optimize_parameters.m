% parameters = [1.03256446674517 0.836262987259303 0.0967320484791088 0.974348858598414];  % tau_1 = .001 beta = .1
% parameters = [2.23731408487711 0.757726476240968 1.09385691270112 104.376090499813];  % tau_1 = .01 beta = .1
% parameters = [4.09653432975792 1.90326374670256 0.873909949641518 484.498224716212]; % tau_1 = .1 beta = .1
% parameters = [42 1000 .00003 .08];
% parameters = [45 1500 .000035 .09];
% parameters = [50 1500 .00002 .101];
% parameters = [55 1700 .000015 .102];
parameters = [.7 5 .2 .1];
options = optimset('Display','iter','PlotFcns',@optimplotfval);
fval = 1;

get_relaxation_functions_difference(parameters);
% while fval > .05
%   [parameters,fval] = fminsearch(@get_relaxation_functions_difference,parameters,options);
% end


%%
function [relaxation_functions_difference] = get_relaxation_functions_difference(parameters)

% stat. model parameters
E0 = parameters(1); % 42 45 50 55
E1 = parameters(2); % 1000 1500 1500 1700
theta = parameters(3); % .00003 .000035 .00002 .000015
tau_1 = parameters(4); % 1e0 1e0
Lo = 1 / (tau_1 * E1);
lambda = 1;
gamma_stat = 1; % /(E1)TODO: неправильно указана зависимость (от E1)

% time
time_step = 0.002442598925256;
time_length = 10;
n_time_points = 1200;
time = linspace(time_step,time_length/2,n_time_points)';
% time_length = 10;
% n_time_points = 2^11;
% time = linspace(0, time_length/2, n_time_points)';
% time_step = time(2) - time(1);

% deformation|identation
eps = .15; % .08 .09 .101 .102
deps_dt = eps / length(time) / time_step;

% the KV-fractional model
E = 1;
tau = 1e-2; % 1e-1 1e-1 1e-1 1e-1
alpha = .75;
beta = .1; % .0 .1 .2 .3
tip_geometry = 2;
F = E * gamma(tip_geometry + 1) * (time / time(end)).^tip_geometry ...
  .* ( 1 / gamma(tip_geometry + 1 - alpha) * (time / tau).^-alpha ...
      + 1 / gamma(tip_geometry + 1 - beta) * (time / tau).^-beta );
relaxation_function_frac_model = E / gamma(1 - alpha) * (time / tau).^-alpha ...
  + E / gamma(1 - beta) * (time / tau).^-beta;
relaxation_function_frac_model = (relaxation_function_frac_model - relaxation_function_frac_model(end)) / (relaxation_function_frac_model(1) - relaxation_function_frac_model(end));
F = [flip(F); zeros(size(F))];

% stat. model
[~,fit_model] = get_xi_approximation;

sigma_0 = 0;
initial_conditions = sigma_0;
% eps_o_0 = (eps * time(1)) - (sigma_0 - E0 * eps * time(1)) / E1;
% initial_conditions = eps_o_0;

% [~,eps_o] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps,deps_dt,E1,0,Lo,0,theta,lambda,gamma_stat,fit_model),time,initial_conditions);
% sigma_o =  E1 * (eps * time - eps_o);
% 
% sigma = E0 * eps * time + sigma_o;

[~,sigma] = ode15s(@(t,sigma) get_rhs_expr(t,sigma,eps,deps_dt,E0,E1,Lo,0,theta,lambda,gamma_stat,fit_model),time,initial_conditions);
hold on;plot(linspace(0,eps,length(sigma)),sigma)
sigma = [flip(sigma); zeros(size(sigma))];

% time
time = [time; linspace(time(end) + time_step,time_length,n_time_points)'];

% wavelet-based decomposition stress-strain curve
n_scales = 2^10;
ssc_scale_size = 10*time_step;
ssc_min_scale = 1;
ssc_max_scale = 1e1;
ssc_scales = exp(log( (ssc_max_scale) / (ssc_min_scale)) * ((1:n_scales)' - 1) / n_scales );

eps_scales = ssc_scales * ssc_scale_size;
[wt_coefficients_0_order,wt_coefficients_1_order,ssc_wt_coefficients_2_order] = get_wavelet_transform_coefficients(sigma,time,ssc_scale_size,n_scales,ssc_max_scale,ssc_min_scale);
mask = ones(size(wt_coefficients_1_order));
mask(1000:end,:) = 0;
wt_coefficients_1_order = abs(wt_coefficients_1_order) .* mask;


[~,ssc_indices] =  max(wt_coefficients_1_order,[],1);
for n_scale = 1:n_scales
  ssc_maxima_line(n_scale) = wt_coefficients_1_order(ssc_indices(n_scale),n_scale);
end
% [ssc_maxima_line,ssc_indices] = max(wt_coefficients_1_order,[],1);
% n_edge_effects_points = 1100;
% mask = ones(size(wt_coefficients_1_order));
% mask(1:n_edge_effects_points,:) = 0;
% ssc_maxima_lines = islocalmax(wt_coefficients_1_order,1) .* mask;
% ssc_indices = zeros(n_scales,1);
% ssc_maxima_line = zeros(n_scales,1);
% for n_scale = 1:n_scales
%   ssc_indices(n_scale) = find(ssc_maxima_lines(:,n_scale) == 1,1);
%   ssc_maxima_line(n_scale) = wt_coefficients_1_order(ssc_indices(n_scale),n_scale);
% end

% wavelet-based decomposition forces-indentation curve
% n_scales = 2^10;
fic_scale_size = 10 * time_step;
fic_min_scale = 1;
fic_max_scale = 1e1;
fic_scales = exp(log( (fic_max_scale) / (fic_min_scale)) * ((1:n_scales)' - 1) / n_scales );
time_scales = fic_scales * fic_scale_size;

[~,~,wt_coefficients_2_order] = get_wavelet_transform_coefficients((F),time,fic_scale_size,n_scales,fic_max_scale,fic_min_scale);
wt_coefficients_2_order(wt_coefficients_2_order < 0) = 0;

[fic_maxima_line,fic_indices] = max(wt_coefficients_2_order,[],1);
% n_edge_effects_points = 790;
% mask = ones(size(wt_coefficients_1_order));
% mask(1:n_edge_effects_points,:) = 0;
% maxima_lines = islocalmax(wt_coefficients_1_order,1) .* mask;
% indices = zeros(n_scales,1);
% maxima_line = zeros(n_scales,1);
% for n_scale = 1:n_scales
%   indices(n_scale) = find(maxima_lines(:,n_scale) == 1,1);
%   maxima_line(n_scale) = wt_coefficients_1_order(indices(n_scale),n_scale);
% end

relaxation_functions_difference = max(abs(fic_maxima_line / max(fic_maxima_line) - ssc_maxima_line / max(ssc_maxima_line)));

% visualization
%
figure(1);
subplot(2,1,1);hold on;
plot(time,sigma,'k','LineWidth',3);
% plot(time,wt_coefficients_0_order(:,1))
ylabel('{\it\sigma}, МПа')
xlabel('{\itt}, с')
set_figure;
subplot(2,1,2);
plot(time,F,'k','LineWidth',3);
ylabel('{\itF}, МПа')
xlabel('{\itt}, с')
set_figure;

figure(2);
% subplot(2,1,1);
% mesh(wt_coefficients_1_order(:,1:end)'); % n_edge_effects_points:end-n_edge_effects_points
imagesc(wt_coefficients_1_order(1:end,1:end)');
% imagesc(abs(wt_coefficients_1_order(1:end,1:end)'));
% mesh(log10(abs(wt_coefficients_1_order(:,1:end)')));
colormap('jet')
axis xy;
hold on; plot(ssc_indices,1:n_scales,'w')
% subplot(2,1,2);
% % mesh(wt_coefficients_2_order(:,1:end)'); % n_edge_effects_points:end-n_edge_effects_points
% imagesc(wt_coefficients_2_order(1:end,1:end)');
% % mesh(log10(abs(wt_coefficients_2_order(:,1:end)')));
% colormap('jet')
% axis xy;
% hold on; plot(fic_indices,1:n_scales,'r')
xlim([0 1000])

% figure(3);
% plot(wt_coefficients_1_order(:,1));

figure(4);hold on;
plot(time_scales,ssc_maxima_line / max(ssc_maxima_line),'k','LineWidth',3); %
% yyaxis right;
plot(time_scales,fic_maxima_line / max(fic_maxima_line),'k:','LineWidth',3); %
ylabel('{\itR}, МПа')
xlabel('{\itt}, с')
legend({'{\itR} при \beta = 0,1','{\itR}^{stat} при \beta = 0,1'})
set_figure
%}
end

function [rhs_expr] = get_rhs_expr(t,varibles,eps,deps_dt,E1,E2,Lo,nu,theta,lambda,gamma,fit_model)
% Функция для расчета правой части дифференциального уравнения

if nargin < 12
  [~,fit_model] = get_xi_approximation;
end
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

% if E2 == 0
  eps_t = eps * t;
  eps_o = eps_t - ( varibles(1) - E1 * eps_t ) / E2;
  ksi = feval(fit_model,eps_o);
  dF_deps_o = 2/3 * theta / gamma * ksi - lambda / gamma * eps_o;
  nu = 1 / Lo;

  rhs_expr(1) = E1 * E2 / nu * eps_t...
  + (E1 + E2) * deps_dt ...
  - E2 / nu * varibles(1) ...
  + E2 / nu * dF_deps_o;
% else
%   rhs_expr(1) = E2 * deps_dt - E2/nu * varibles(1);
% end

end


