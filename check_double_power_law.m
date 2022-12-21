
% set number of force curve and spans of model parameters
n_force_curve = 4; % 1..6
E1_span = [15; 90; 450; 700; 900; 1700];
tau1_span = [1; 1; 1; 2; 3; 2];
theta_span = [.16; .16; .16; .16; .16; .16];
eps_loading_span = [.2; .22; .25; .3; .33; .39];

% time
time_step = .0013;
time = (time_step:time_step:10)';
time_loading = time(time <= 1);
time_dwell = time(time > 1);

% deformation|identation
eps_loading = eps_loading_span(n_force_curve);
eps_dwell = .1;
deps_dt = eps_loading * time_step;
deps_dt_dwell = 0;

% model parameters
E1 = E1_span(n_force_curve);
tau_1 = tau1_span(n_force_curve);
Lo = 1 / (tau_1 * E1);
lambda = 1;
gamma = 1/(E1);

% initial conditions
F = readmatrix('data.txt') * 1000;
F_loading = F(time <= 1,:);
F_dwell = F(time > 1,:);
sigma_0_loading = F_loading(1,n_force_curve);
sigma_0_dwell = abs(F_dwell(1,n_force_curve) - F_dwell(end,n_force_curve));
eps_o_0_loading = 0 - sigma_0_loading / E1;
eps_o_0_dwell = eps_dwell - sigma_0_dwell / E1;

% model parameters
E0_loading = F_dwell(end,n_force_curve) / eps_dwell; % F_loading(end,n_force_curve) / eps_loading
E0_dwell = F_dwell(end,n_force_curve) / eps_dwell;

% effective temperature factor
theta = theta_span(n_force_curve);

% loading regime
initial_conditions = 0;
[~,eps_o_loading] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps_loading,deps_dt,E1,0,Lo,0,theta,lambda,gamma),time_loading,initial_conditions);
sigma_o_loading =  E1 * (eps_loading * time_loading - eps_o_loading);

% dwell regime
initial_conditions = eps_o_0_dwell;
[~,eps_o_dwell] = ode15s(@(t,eps_o) get_rhs_expr(t,eps_o,eps_dwell,deps_dt_dwell,E1,0,Lo,0,theta,lambda,gamma),time_dwell,initial_conditions);
sigma_o_dwell =  E1 * (eps_dwell - eps_o_dwell);

% visualization results
figure(1);hold on
% plot(time_loading,sigma_o_loading,'r','LineWidth',2)
plot(...
  log10(time_loading),...
  log10((E0_loading * eps_loading * time_loading + sigma_o_loading) ./ (eps_loading * time_loading) .* (eps_loading * time_loading).^(3/2)),... % 
  'k','LineWidth',2)
plot(log10(time_loading),log10(F_loading(:,n_force_curve)),'k--','LineWidth',2)
xlabel('time');
ylabel('\sigma');
set_figure;

figure(2);hold on
%   plot(time_dwell,sigma_o_dwell,'r','LineWidth',2)
plot((time_dwell),(E0_dwell * eps_dwell * ones(size(time_dwell)) + sigma_o_dwell),'k','LineWidth',2)
plot((time_dwell),(F_dwell(:,n_force_curve)),'k--','LineWidth',2)
xlabel('time');
ylabel('\sigma');
set_figure;

