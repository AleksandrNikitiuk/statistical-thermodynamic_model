%% Оценка периода собственных колебаний сегмента и пачки актиновых филаментов
m = 1.66053906660e-24;                                    % в кг
% e = 1.602176634e-19;                                      % в Дж
e = 1e0;                                                  % в Дж 
l = 1e0;                                                 % в м
tau_actin_segment = sqrt((m*10) * (l*10)^2 / (e*10))      % в сек
tau_actin_bundles = sqrt((m*100) * (l*1000)^2 / (e*100))  % в сек

%% Решение уравнения самосогласования
effective_field_shear_case = (-10:.11:30)';
effective_field_orientation_case = (-30:.011:30)';

is_dimensionless = 1;
[shear_deformation] = get_shear_deformation(effective_field_shear_case,[],is_dimensionless);
[orientation_deformation] = get_orientation_deformation(effective_field_orientation_case,[],is_dimensionless);

% Аппроксимация эффективных полей
fit_type_orientation_case = fittype('rat34');
coefficients_orientation_case = [-1.142e+05 -8815; 7.625e+04 -2.067e+04; 1.232e+05 1.677e+04; 255.3 1.24e+04];
[effective_field_orientation_case_approximation,fit_model_orientation_case] = ...
  get_effective_field_approximation(...
  effective_field_orientation_case,...
  orientation_deformation,...
  fit_type_orientation_case,...
  coefficients_orientation_case);

effective_field_shear_case = ...
  effective_field_shear_case(shear_deformation > 0) - shear_deformation(shear_deformation > 0);
shear_deformation = shear_deformation(shear_deformation > 0);
fit_type_shear_case = fittype('rat34');
coefficients_shear_case = [16.77 -3.353; -83.97 6.894; 327.2 34.51; -1.489 34.54];
[effective_field_shear_case_approximation,fit_model_shear_case] = ...
  get_effective_field_approximation(...
  effective_field_shear_case,...
  shear_deformation,...
  fit_type_shear_case,...  % 
  coefficients_shear_case); % 

figure(1);hold on;
plot(shear_deformation,effective_field_shear_case,'k');
plot(shear_deformation,effective_field_shear_case_approximation,'r--',...
  'LineWidth',3);
xlabel(['{\it' char(949) '}_{\its}' ]);ylabel('{\it\xi_s}')

figure(2);hold on;
plot(orientation_deformation,effective_field_orientation_case,'k');
plot(orientation_deformation,effective_field_orientation_case_approximation,'r--',...
  'LineWidth',3);
xlabel(['{\it' char(949) '}_{\ito}' ]);ylabel('{\it\xi_o}')

%% Свободная энергия ориентационная часть
orientation_deformation = -.3:.01:.8;
sigmas = [.1 .2 .3]; % -.7:.01:.4
chi_orientation = .1;

effective_field_integral_orientation_case = ...
  get_effective_field_integral_orientation_case(orientation_deformation);

F_orientation = zeros(numel(effective_field_integral_orientation_case),numel(sigmas));
for sigma = sigmas
  F_orientation(:,sigma == sigmas) = ...
    effective_field_integral_orientation_case ...
    - sigmas(sigma == sigmas) .* orientation_deformation ...
    - chi_orientation^-1 * orientation_deformation.^2 / 2;
  F_orientation(:,sigma == sigmas) = ...
    F_orientation(:,sigma == sigmas) - mean(F_orientation(:,sigma == sigmas));

  F_orientation_barriers = get_energy_barrier(F_orientation(:,sigma == sigmas))

  figure(3);hold on;
  plot(orientation_deformation,F_orientation(:,sigma == sigmas),... %  - max(F_orientation(:,sigma == sigmas))
    'LineWidth',3);
  xlabel(['{\it' char(949) '}_{\ito}' ]);ylabel('\Psi_{\ito} - <\Psi_{\ito}>')
end
axis tight;
colororder('k'); ax = gca; ax.LineStyleOrder = ["--";"-";  ":"];

%% Свободная энергия сдвиговая часть
shear_deformation = 0:.01:25;
sigma = 3;
chi_shear = [1.25 1.16 1.07];

effective_field_integral_shear_case = ...
  get_effective_field_integral_shear_case(shear_deformation);

F_shear = zeros(numel(effective_field_integral_shear_case),numel(chi_shear));
for chi = chi_shear
  F_shear(:,chi==chi_shear) = ...
    effective_field_integral_shear_case  + shear_deformation.^2 / 2 ...
    - sigma .* shear_deformation - chi_shear(chi==chi_shear)^-1 * shear_deformation.^2 / 2;
  F_shear(:,chi==chi_shear) = ...
    F_shear(:,chi==chi_shear) - mean(F_shear(:,chi==chi_shear));
  
  F_shear_barriers = get_energy_barrier(F_shear(:,chi==chi_shear))

  figure(4);hold on;
  plot(shear_deformation,F_shear(:,chi==chi_shear),'LineWidth',3);
  xlabel(['{\it' char(949) '}_{\its}' ]);ylabel('\Psi_{\its}')
end
axis tight;
% colororder('k'); ax = gca; ax.LineStyleOrder = ["-"; "--"; ":"];

%% Характерные зависимости параметра порядка, описывающего ориентирование сегментов филаментов, от приложенного напряжения
chi_orientation = .2; % .06 .1 .2
orientation_deformation = (-.3:.01:.8)';
sigma_orientation = (feval(fit_model_orientation_case,orientation_deformation) - chi_orientation^-1 * orientation_deformation);

figure(1);hold on;
plot(sigma_orientation,orientation_deformation,'LineWidth',3);
xlabel('{\it\sigma}');
ylabel(['{\it' char(949) '_o}']);

%% Характерные зависимости параметра порядка, описывающего скольжение пачек филаментов, от приложенного напряжения
shear_deformation = (0:.01:25)';
chi_shear = .9; % [1.25 1.16 1.07 1.01 .9]
sigma_shear = (...
  feval(fit_model_shear_case,shear_deformation) ...
  + shear_deformation ...
  - chi_shear.^-1 * shear_deformation);

figure(1);hold on;
plot(sigma_shear,shear_deformation,'LineWidth',3);
xlabel('{\it\sigma}');
ylabel(['{\it' char(949) '_s}']);


%%

% Time
n_time_steps = 8e4;
time_step = 1e-4;                                                           % с

% Varibles initialization
sigma_shear = zeros(n_time_steps,1);
sigma_reversible = zeros(size(sigma_shear));
sigma_elastic = zeros(size(sigma_shear));
orientation_strain_increment = zeros(size(sigma_shear));
chi_shear_increment = zeros(size(sigma_shear));
shear_strain_increment = zeros(size(sigma_shear));
strain = zeros(size(sigma_shear));
orientation_strain = zeros(size(sigma_shear));
chi_shear = zeros(size(sigma_shear));
shear_strain = zeros(size(sigma_shear));
orientation_strain_rate = zeros(size(sigma_shear));
shear_strain_rate = zeros(size(sigma_shear));
elastic_strain = zeros(size(sigma_shear));
sigma_elastic_increment = zeros(size(sigma_shear));
sigma_reversible_increment = zeros(size(sigma_shear));
strain_rate = zeros(size(sigma_shear));
orientation_relaxation_time = zeros(size(sigma_shear));
shear_relaxation_time = zeros(size(sigma_shear));

% Initial conditions
sigma_shear(:) = 1.0e0;                                                       % Па; 1 и 6 Па - демонстрация ориентационной и микросдвиговой деформации соответственно
chi_shear(1) = 1.3;                                                        % 1 [1.3 1.1 1.01]

% Model parameters for orientation case
orientation_deformation = -.3:.01:.8;                                       % 1
effective_field_integral_orientation_case = ...
    get_effective_field_integral_orientation_case(orientation_deformation);
F_orientation = zeros(n_time_steps,numel(orientation_deformation));
k_boltzmann = 1.380649e-23;                                                 % Дж/К
temperature = 677.5;                                                          % К
theta = k_boltzmann * temperature / 2;
lambda_orient = 2e1;                                                        % Па
gamma_orient = 1e-21;                                                       % м^3
chi_orientation = theta / (lambda_orient * gamma_orient);
initial_orientation_relaxation_time = .5086e-2;                            % c

% Model parameters for shear case
shear_deformation = 0:.01:25;                                               % 1
effective_field_derivative_shear_case = ...
    get_effective_field_integral_shear_case(shear_deformation);
F_shear = zeros(n_time_steps,numel(shear_deformation));
initial_shear_relaxation_time = 1.2886e-1;                                  % c
alpha = 1e21;                                                               % Н / м^5
bundles_concentration = 1e19;                                               % м^-3

% Other model parameters
chi_shear_kin_coef = 1e-3;                                                   % ?
vortex = zeros(3); vortex(1,2) = 1; vortex(2,1) = 1;
I = eye(3);
[C_I,C_II] = get_fourth_rank_isotropic_tensors;
cell_shear_modulus = 1e2;                                                   % Па
actin_shear_modulus = lambda_orient;                                        % Па
C1 = actin_shear_modulus / 2;                                               % Па
poissons_ratio = .49;
D1 = 0; % 2 * actin_shear_modulus * poissons_ratio / (1 - 2 * poissons_ratio)
C1_reversible = cell_shear_modulus / 2;                                     % Па
D1_reversible = 0; % 2 * cell_shear_modulus * poissons_ratio / (1 - 2 * poissons_ratio)

% Initial kinematics
F = I;
B = zeros(n_time_steps,3,3);
B(1,:,:) = tensorprod(F,F',2,1);
B(:,3,3) = 1;
B_elastic = zeros(n_time_steps,3,3);
B_elastic(:,3,3) = 1;
B_elastic(1,:,:) = tensorprod(F,F',2,1);

% Approximations of effective fileds
effective_field_shear_case = (-10:.11:30)';
effective_field_orientation_case = (-30:.011:30)';

is_dimensionless = 1;
[s_deformation] = get_shear_deformation(effective_field_shear_case,[],is_dimensionless);
[o_deformation] = get_orientation_deformation(effective_field_orientation_case,[],is_dimensionless);

fit_type_orientation_case = fittype('rat34');
coefficients_orientation_case = [-1.142e+05 -8815; 7.625e+04 -2.067e+04; 1.232e+05 1.677e+04; 255.3 1.24e+04];
[~,fit_model_orientation_case] = ...
  get_effective_field_approximation(...
  effective_field_orientation_case,...
  o_deformation,...
  fit_type_orientation_case,...
  coefficients_orientation_case);

effective_field_shear_case = ...
  effective_field_shear_case(s_deformation > 0) - s_deformation(s_deformation > 0);
s_deformation = s_deformation(s_deformation > 0);
fit_type_shear_case = fittype('rat34');
coefficients_shear_case = [16.77 -3.353; -83.97 6.894; 327.2 34.51; -1.489 34.54];
[~,fit_model_shear_case] = ...
  get_effective_field_approximation(...
  effective_field_shear_case,...
  s_deformation,...
  fit_type_shear_case,...  % 
  coefficients_shear_case); %

for t = 1:(n_time_steps - 1)
  
  % get orient_kin_coef and shear_kin_coef values
  F_orientation(t,:) = ...
    effective_field_integral_orientation_case ...
    - gamma_orient / theta * sigma_elastic(t) .* orientation_deformation ...
    - chi_orientation^-1 * orientation_deformation.^2 / 2;
  F_orientation(t,:) = F_orientation(t,:) - mean(F_orientation(t,:));
  F_orientation_barrier = get_energy_barrier(F_orientation(t,:));
  orientation_relaxation_time(t) = initial_orientation_relaxation_time * ...
    exp(F_orientation_barrier);
  orient_kin_coef = 1 / (3 * actin_shear_modulus * orientation_relaxation_time(t)); % 
  
  F_shear(t,:) = ...
    effective_field_derivative_shear_case  + shear_deformation.^2 / 2 ...
    - gamma_orient / theta * sigma_elastic(t) .* shear_deformation ...
    - chi_shear(t)^-1 * shear_deformation.^2 / 2;
  F_shear(t,:) = F_shear(t,:) - mean(F_shear(t,:));
  F_shear_barrier = get_energy_barrier(F_shear(t,:));
  shear_relaxation_time(t) = initial_shear_relaxation_time ...
    * exp((F_shear_barrier)); % F_orientation_barrier + ?
  shear_kin_coef = 1 / (3 * cell_shear_modulus * shear_relaxation_time(t)); %  ? shear_modulus
  % (end) get orient_kin_coef and shear_kin_coef values]

  % get derivatives
  dimensionless_orientation_strain = orientation_strain(t);
  dF_orientation_strain = theta / gamma_orient * (...
    feval(fit_model_orientation_case,dimensionless_orientation_strain)...
    - chi_orientation^-1 * dimensionless_orientation_strain...
    );
  
  dimensionless_shear_deformation = sqrt(alpha / theta) / bundles_concentration * shear_strain(t);
  dF_shear_deformation = sqrt(alpha * theta) * (...
    feval(fit_model_shear_case,dimensionless_shear_deformation)...
    + dimensionless_shear_deformation...
    - chi_shear(t).^-1 * dimensionless_shear_deformation...
    );

  dF_chi_shear = (2*chi_shear(t).^2).^-1 .* shear_strain(t).^2; % heaviside(orientation_strain(t) - .5) * 
  % end get derivatives

  orientation_strain_increment(t) = ...
    orient_kin_coef * (sigma_elastic(t) - dF_orientation_strain) * time_step;
  chi_shear_increment(t) = - chi_shear_kin_coef * dF_chi_shear * time_step;
  shear_strain_increment(t) = ...
    shear_kin_coef * (sigma_elastic(t) - dF_shear_deformation) * time_step;

  orientation_strain_rate(t) = orientation_strain_increment(t) / time_step;
  shear_strain_rate(t) = shear_strain_increment(t) / time_step;

  strain_rate(t) = (...
    B_elastic(t,1,2) / time_step + B(t,1,2) / time_step ...
    - sigma_shear(t + 1) / (2 * C1 * time_step) ...
    + 2 * (orientation_strain_rate(t) + shear_strain_rate(t)) ...
    * (B_elastic(t,2,2) + B_elastic(t,1,2)) ...
    ) / ...
    (2 * (B(t,1,2) + B(t,2,2) + B_elastic(t,1,2) + B_elastic(t,2,2)));

  vortex = strain_rate(t) / 2 ...
    * ([1 0 0]' .* [0 1 0] - [0 1 0]' .* [1 0 0]);
  B(t + 1,1,1) = B(t,1,1) + time_step * (...
    -vortex(1,2) * ( B(t,2,1) + B(t,1,2) ) - strain_rate(t) * ( B(t,1,2) + B(t,2,1) ));
  B(t + 1,2,2) = B(t,2,2) + time_step * (...
    -vortex(1,2) * ( B(t,2,1) + B(t,1,2) ) - strain_rate(t) * ( B(t,1,2) + B(t,2,1) ));
  B(t + 1,1,2) = B(t,1,2) + time_step * (...
    vortex(1,2) * ( B(t,1,1) - B(t,2,2) ) - strain_rate(t) * ( B(t,1,1) + B(t,2,2) ));
  B(t + 1,2,1) = B(t,2,1) + time_step * (...
    vortex(2,1) * ( B(t,2,2) - B(t,1,1) ) - strain_rate(t) * ( B(t,2,2) + B(t,1,1) ));

  vortex_elastic = (strain_rate(t) - orientation_strain_rate(t) - shear_strain_rate(t)) / 2 ...
    * ([1 0 0]' .* [0 1 0] - [0 1 0]' .* [1 0 0]);
  B_elastic(t + 1,1,1) = B_elastic(t,1,1) + time_step * (...
    -vortex_elastic(1,2) * ( B_elastic(t,2,1) + B_elastic(t,1,2) ) - (strain_rate(t) - orientation_strain_rate(t) - shear_strain_rate(t)) * ( B_elastic(t,1,2) + B_elastic(t,2,1) ));
  B_elastic(t + 1,2,2) = B_elastic(t,2,2) + time_step * (...
    -vortex_elastic(2,1) * ( B_elastic(t,1,2) + B_elastic(t,2,1) ) - (strain_rate(t) - orientation_strain_rate(t) - shear_strain_rate(t)) * ( B_elastic(t,2,1) + B_elastic(t,1,2) ));
  B_elastic(t + 1,1,2) = B_elastic(t,1,2) + time_step * (...
    vortex_elastic(1,2) * ( B_elastic(t,1,1) - B_elastic(t,2,2) ) - (strain_rate(t) - orientation_strain_rate(t) - shear_strain_rate(t)) * ( B_elastic(t,1,1) + B_elastic(t,2,2) ));
  B_elastic(t + 1,2,1) = B_elastic(t,2,1) + time_step * (...
    vortex_elastic(2,1) * ( B_elastic(t,2,2) - B_elastic(t,1,1) ) - (strain_rate(t) - orientation_strain_rate(t) - shear_strain_rate(t)) * ( B_elastic(t,2,2) + B_elastic(t,1,1) ));
  
  sigma_elastic(t + 1) = 2 * C1 * B_elastic(t + 1,1,2);
  sigma_reversible(t + 1) = 2 * C1 * B(t + 1,1,2);

  orientation_strain(t + 1) = orientation_strain(t) + orientation_strain_increment(t);
  chi_shear(t + 1) = chi_shear(t) + chi_shear_increment(t);
  shear_strain(t + 1) = shear_strain(t) + shear_strain_increment(t);

end

eps = 1e-1;
time = time_step * (1:n_time_steps);
% time = time(abs(diff(strain_rate)) < eps);
% orientation_strain_rate = orientation_strain_rate(abs(diff(strain_rate)) < eps);
% shear_strain_rate = shear_strain_rate(abs(diff(strain_rate)) < eps);
% chi_shear = chi_shear(abs(diff(strain_rate)) < eps);
% strain_rate = strain_rate(abs(diff(strain_rate)) < eps);
strain = cumtrapz(time,strain_rate);

mean_square_error = ...
  mean( sqrt( ...
  abs(sigma_elastic(2:end) + sigma_reversible(2:end) - sigma_shear(2:end)).^2 ...
  ./ abs(sigma_shear(2:end)).^2 ) );


line_width = 2;

figure(7);hold on;
plot(orientation_deformation,F_orientation(1,:),'b','LineWidth',line_width)
plot(orientation_deformation(31),F_orientation(1,31),'bo','MarkerSize',7,'LineWidth',line_width)
plot(orientation_deformation,F_orientation(100,:),'r','LineWidth',line_width)
plot(orientation_deformation(52),F_orientation(100,52),'ro','MarkerSize',7,'LineWidth',line_width)
xlabel(['{\it' char(949) '_o}']);
ylabel('\Psi_{\ito} - <\Psi_{\ito}>');
legend({ ...
  '\Psi_{\ito}(0.003,{\itε_o}({\itt}))', ...
  '\Psi_{\ito}(0.003,{\itε_o}(0.003))', ...
  '\Psi_{\ito}(0.005,{\itε_o}({\itt}))', ...
  '\Psi_{\ito}(0.005,{\itε_o}(0.005))'},'Location','best');
 
figure(6);hold on;
plot(shear_deformation,F_shear(1100,:),'b','LineWidth',line_width)
plot(shear_deformation(1),F_shear(1100,1),'bo','MarkerSize',7,'LineWidth',line_width)
plot(shear_deformation,F_shear(1200,:),'r','LineWidth',line_width)
plot(shear_deformation(50),F_shear(1200,50),'ro','MarkerSize',7,'LineWidth',line_width)
plot(shear_deformation,F_shear(1250,:),'m','LineWidth',line_width)
plot(shear_deformation(100),F_shear(1250,100),'mo','MarkerSize',7,'LineWidth',line_width)
xlabel(['{\it' char(949) '_s}']);
ylabel('\Psi_{\its} - <\Psi_{\its}>')
legend({ ...
  '\Psi_{\its}(0.11,{\itε_s}({\itt}))', ...
  '\Psi_{\its}(0.11,{\itε_s}(0.11))', ...
  '\Psi_{\its}(0.12,{\itε_s}({\itt}))', ...
  '\Psi_{\its}(0.12,{\itε_s}(0.12))', ...
  '\Psi_{\its}(0.125,{\itε_s}({\itt}))', ...
  '\Psi_{\its}(0.125,{\itε_s}(0.125))'}, ...
  'Location','best')

figure(5);hold on;
colororder({'b','r'})
yyaxis left;
plot(time(2:end-1),orientation_relaxation_time(2:(end-1)),'LineWidth',line_width)
ylabel('{\it\tau_o}, s')
yyaxis right;
plot(time(2:end-1),shear_relaxation_time(2:(end-1)),'LineWidth',line_width)
xlabel('{\itt}, s');
ylabel('{\it\tau_s}, s')

figure(4);hold on;
plot(time(2:end-1),strain_rate(2:end-1),'k','LineWidth',line_width + 1)
plot(time(2:end-1),orientation_strain_rate(2:end-1),'b','LineWidth',line_width)
plot(time(2:end-1),shear_strain_rate(2:end-1),'r','LineWidth',line_width)
plot(time(2:end-1),strain_rate(2:end-1) - orientation_strain_rate(2:end-1) - shear_strain_rate(2:end-1),'g','LineWidth',line_width)
plot(time(2:end-1),strain_rate(2:end-1) - orientation_strain_rate(2:end-1) - shear_strain_rate(2:end-1)  + orientation_strain_rate(2:end-1) + shear_strain_rate(2:end-1),'m--','LineWidth',line_width)
xlabel('{\itt}, s');
ylabel('{\itD}, s^{-1}');
legend({'{\itD}','{\itD_o}','{\itD_s}','{\itD_e}','{\itD_e} + {\itD_o} + {\itD_s}'},'Location','best'); % 

figure(3);hold on;
plot((2:n_time_steps)*time_step,sigma_shear(2:end),'k','LineWidth',line_width + 1)
plot((2:n_time_steps)*time_step,sigma_reversible(2:end),'b','LineWidth',line_width)
plot((2:n_time_steps)*time_step,sigma_elastic(2:end),'r','LineWidth',line_width)
plot((2:n_time_steps)*time_step,sigma_reversible(2:end) + sigma_elastic(2:end),'m--','LineWidth',line_width)
xlabel('{\itt}, s');
ylabel('{\it\sigma}, Pa')
legend({'{\it\sigma_s}','{\it\sigma_r}','{\it\sigma_e}','{\it\sigma_r} + {\it\sigma_e}'},'Location','best');

figure(2);hold on;
plot((1:n_time_steps)*time_step,chi_shear)
xlabel('{\itt}, s');
ylabel('{\it\chi_s}');

figure(1);hold on;
plot(time,strain + abs(min(strain)),'k','LineWidth',line_width)
plot(time,orientation_strain,'b','LineWidth',line_width)
plot(time,shear_strain,'r','LineWidth',line_width)
xlabel('{\itt}, s');
ylabel(['{\it' char(949) '}'])
legend({['{\it' char(949) '}'],['{\it' char(949) '_o}'],['{\it' char(949) '_s}']},'Location','best');

%% Кинетика актиновых филаментов, обусловленная ориентированием их сегментов
sigma = 5;
chi_orientation = .3;
G = 1;
tau = 1;

time = linspace(0,1,2^10);
initial_conditions = 0;

effective_field_orientation_case = (-10:.011:10)';
[orientation_deformation] = get_orientation_deformation(effective_field_orientation_case,[],1);
fit_type_orientation_case = fittype('rat34');
coefficients_orientation_case = [-6.737e+04 6454; -1.31e+04 -1.429e+05; 4.14e+05 -1.466e+04; -2317 3.194e+05];
[~,fit_model_orientation_case] = ...
  get_effective_field_approximation(...
  effective_field_orientation_case,...
  orientation_deformation,...
  fit_type_orientation_case,...
  coefficients_orientation_case);

[~,orientation_deformation] = ode15s(@(t,orientation_deformation) ...
  get_orientation_deformation_increment(t,orientation_deformation,fit_model_orientation_case,sigma,chi_orientation,G,tau),...
  time,initial_conditions);

deformation = -.49:.01:.99;
effective_field_integral_orientation_case = ...
  get_effective_field_integral_orientation_case(deformation);

F_orientation = zeros(numel(effective_field_integral_orientation_case),numel(time));
F_orientation_barriers = zeros(numel(time),1);
for t = time
  F_orientation(:,t == time) = ...
    effective_field_integral_orientation_case ...
    - sigma * t .* deformation ...
    - chi_orientation^-1 * deformation.^2;

  F_orientation_barriers(t == time) = get_energy_barrier(F_orientation(:,t == time));

end

check = 0;
figure('Color','w'); % ,'units','normalized','outerposition',[0 0 1 1]
for k = 1:numel(time)  
  plot(deformation, F_orientation(:,k),'k-');
  [~,i] = min(abs(deformation - orientation_deformation(k)));
  hold on;
  plot(deformation(i),F_orientation(i,k),'ko')
  hold off;
  xlabel(['{\it' char(949) '_o}']);
  ylabel('{\itF_o}');
  axis tight;
%   ylim([-10^-0 * A 10^-0 *A]);
%   title(['t = ' num2str(T(k)) ' t.u.']);
  f = getframe(gcf);
  if (check == 0)
    [im,map] = rgb2ind(f.cdata,2^10,'nodither');
    check = 1;
  else
    im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  end  
end
% 
% % gif_file_name = 'breathers_T_nonzero_1_dynamics.gif'; % on-site_breather inter-site_breather breathers
% % imwrite(im,map,gif_file_name,'DelayTime',0,'LoopCount',inf);w(2);

% figure(3);hold on;
% plot(sigma*time,orientation_deformation);
% xlabel('\sigma');ylabel(['{\it' char(949) '}']);
% 
% figure(4);hold on
% plot(time,F_orientation_barriers);
% xlabel('{\itt}');ylabel('{\itF_{barrier}}');

%% Кинетика сегментов и пачек актиновых филаментов
sigma = 1e0;                      % Па
chi_orientation = .3;             % 1
chi_shear = 1.04;                 % 1
G = 1e2;                          % Па
tau_orientation = 1.2886e-11;     % c
tau_shear = 1.2886e-9;            % c

theta = 1e-21;                    % Дж
alpha = 1e21;                     % Н / м^5
gamma_o = 1e-21;                  % м^3
lamda_o = 1e9;                    % Па
segments_concentration = 1e24;    % м^-3
bundles_concetration = 1e21;      % м^-3 ?

kinetic_coefficient_chi = 1e5;        % ?

time_end = 2e-6;                  % с
time = linspace(0,time_end,2^9);  % c
initial_conditions = [.0; .0; chi_shear];

% sigma_law = @(time,sigma) sigma;
% sigma_law = @(time,sigma) sigma * time;
% sigma_law = @(time,sigma) sigma / ceil(time_end/2) * (time .^ (time < ceil(time_end/2)));
% sigma_law = @(time,sigma) sigma * time.^2;

is_dimensionless = 1;

effective_field_orientation_case = (-10:.011:10)';
[orientation_deformation] = get_orientation_deformation(effective_field_orientation_case,[],is_dimensionless);
fit_type_orientation_case = fittype('rat34');
coefficients_orientation_case = [-6.737e+04 6454; -1.31e+04 -1.429e+05; 4.14e+05 -1.466e+04; -2317 3.194e+05];
[~,fit_model_orientation_case] = ...
  get_effective_field_approximation(...
  effective_field_orientation_case,...
  orientation_deformation,...
  fit_type_orientation_case,...
  coefficients_orientation_case);

effective_field_shear_case = (-10:.11:30)';
[shear_deformation] = get_shear_deformation(effective_field_shear_case,[],is_dimensionless);
effective_field_shear_case = ...
  effective_field_shear_case(shear_deformation > 0) - shear_deformation(shear_deformation > 0);
shear_deformation = shear_deformation(shear_deformation > 0);
fit_type_shear_case = fittype('rat34');
coefficients_shear_case = [2.024 -2.656; -5.848 4.144; 9.354 -1.004; -0.01012 4.446];
[~,fit_model_shear_case] = ...
  get_effective_field_approximation(...
  effective_field_shear_case,...
  shear_deformation,...
  fit_type_shear_case,...  % 
  coefficients_shear_case); % 

 [~,deformation] = ode15s(@(t,deformation) ...
  get_deformation_increments(...
  t,deformation,fit_model_orientation_case,fit_model_shear_case,sigma,...
  chi_orientation,G,tau_orientation,tau_shear,time_end,theta,alpha,...
  gamma_o,bundles_concetration,kinetic_coefficient_chi),...
  time,initial_conditions);

orientation_deformation = -.99:.01:.99;
effective_field_integral_orientation_case = ...
  get_effective_field_integral_orientation_case(orientation_deformation);
F_orientation = zeros(numel(effective_field_integral_orientation_case),numel(time));
F_orientation_barriers = zeros(numel(time),1);
for t = time
  F_orientation(:,t == time) = ...
    effective_field_integral_orientation_case ...
    - gamma_o / theta * get_sigma(sigma,t,time_end) .* orientation_deformation ...
    - chi_orientation^-1 * orientation_deformation.^2;

  F_orientation_barriers(t == time) = get_energy_barrier(F_orientation(:,t == time));
%   F_orientation_barriers(t == time) = get_current_energy_barrier(F_orientation(:,t == time),orientation_deformation,deformation(t == time,1));

end

shear_deformation = 0:.01:5;
effective_field_integral_shear_case = ...
  get_effective_field_integral_shear_case(shear_deformation);
F_shear = zeros(numel(effective_field_integral_shear_case),numel(time));
F_shear_barriers = zeros(numel(time),1);
for t = time
  F_shear(:,t==time) = ...
    effective_field_integral_shear_case  + shear_deformation.^2 ...
    - get_sigma(sigma,t,time_end) / sqrt(alpha * theta) .* shear_deformation - deformation(t==time,3)^-1 * shear_deformation.^2;
  
  F_shear_barriers(t == time) = get_energy_barrier(F_shear(:,t == time));
%   F_shear_barriers(t == time) = get_current_energy_barrier(F_shear(:,t == time),shear_deformation,deformation(t == time,2));
end
F_shear_barriers(abs(F_shear_barriers) < .1) = 0;

sigma_orientation = (feval(fit_model_orientation_case,orientation_deformation') - chi_orientation^-1 * orientation_deformation');
sigma_shear = zeros(numel(shear_deformation),numel(time));
for t = time
  sigma_shear(:,t == time) = (feval(fit_model_shear_case,shear_deformation') + shear_deformation' - deformation(t == time,3).^-1 * shear_deformation');
end

check = 0;
figure('Color','w'); % ,'units','normalized','outerposition',[0 0 1 1]
pos = get(gcf,'Position');
pos(1:2) = pos(1:2) - 1.5*floor(pos(1:2)/2);pos(3:4) = pos(3:4) * 2;
set(gcf,'Position',pos);
for k = 1:numel(time)  

  subplot(2,2,1);
  plot(orientation_deformation, F_orientation(:,k),'k-');
  [~,i] = min(abs(orientation_deformation - deformation(k,1)));
  hold on;
  plot(orientation_deformation(i),F_orientation(i,k),'ko')
  hold off;
  xlabel(['{\it' char(949) '_o}''']);
  ylabel('{\itF_o}''');
  axis tight;

  subplot(2,2,3);
  plot(sigma_orientation,orientation_deformation,'k');
  sigma_orientation_current = (feval(fit_model_orientation_case,orientation_deformation(i)) - chi_orientation^-1 * orientation_deformation(i));
  hold on;
  plot(sigma_orientation_current,orientation_deformation(i),'ko');
  hold off;
  ylabel(['{\it' char(949) '_o}''']);
  xlabel('{\it\sigma}''');
  axis tight;

  subplot(2,2,2);
  plot(shear_deformation, F_shear(:,k),'k-');
  [~,i] = min(abs(shear_deformation - sqrt(alpha/theta)/bundles_concetration * deformation(k,2)));
  hold on;
  plot(shear_deformation(i),F_shear(i,k),'ko')
  hold off;
  xlabel(['{\it' char(949) '_s}''']);
  ylabel('{\itF_s}''');
  axis tight;

  subplot(2,2,4);
  plot(sigma_shear(:,k),shear_deformation,'k');
  sigma_shear_current = (feval(fit_model_shear_case,shear_deformation(i)) + shear_deformation(i) - deformation(k,3).^-1 * shear_deformation(i));
  hold on;
  plot(sigma_shear_current,shear_deformation(i),'ko');
  hold off;
  ylabel(['{\it' char(949) '_s}''']);
  xlabel('{\it\sigma}''');
  axis tight;

  f = getframe(gcf);
  if (check == 0)
    [im,map] = rgb2ind(f.cdata,2^10,'nodither');
    check = 1;
  else
    im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  end  
end
% gif_file_name = 'breathers_T_nonzero_1_dynamics.gif'; % on-site_breather inter-site_breather breathers
% imwrite(im,map,gif_file_name,'DelayTime',0,'LoopCount',inf);w(2);

% check = 0;
% figure('Color','w'); % ,'units','normalized','outerposition',[0 0 1 1]
% for k = 1:numel(time)  
%   plot(orientation_deformation, F_orientation(:,k),'k-');
%   [~,i] = min(abs(orientation_deformation - deformation(k,1)));
%   hold on;
%   plot(orientation_deformation(i),F_orientation(i,k),'ko')
%   hold off;
%   xlabel(['{\it' char(949) '_o}']);
%   ylabel('{\itF_o}');
%   axis tight;
% %   ylim([-10^-0 * A 10^-0 *A]);
% %   title(['t = ' num2str(T(k)) ' t.u.']);
%   f = getframe(gcf);
%   if (check == 0)
%     [im,map] = rgb2ind(f.cdata,2^10,'nodither');
%     check = 1;
%   else
%     im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
%   end  
% end
% 
% % gif_file_name = 'breathers_T_nonzero_1_dynamics.gif'; % on-site_breather inter-site_breather breathers
% % imwrite(im,map,gif_file_name,'DelayTime',0,'LoopCount',inf);w(2);

sigma_t = zeros(size(time));
for t = time
  sigma_t(t == time) = get_sigma(sigma,t,time_end);
end

figure(3);hold on;

figure(3);hold on;
% plot(time,sigma_t' / (2 * G));
plot(time,deformation(:,1));
plot(time,deformation(:,2));
plot(time,sigma_t' / (2 * G) + deformation(:,1) + deformation(:,2));
xlabel('{\itt}, s');ylabel(['{\it' char(949) '}']);
  
figure(4);hold on;
plot(time,deformation(:,3));
xlabel('{\itt}, s');ylabel('{\it\chi_s}');

figure(5);hold on
plot(time,(G * tau_orientation * exp(F_orientation_barriers)).^-1);
plot(time,(G * tau_shear * exp(F_shear_barriers)).^-1);
% plot(time,(G * 1 * exp(F_orientation_barriers ./ F_shear_barriers)).^-1);
xlabel('{\itt}');ylabel('Kinetic coefficients');%ylabel('{\it\tau}');%ylabel('{\itF_{barrier}}');

%%
function [C_I,C_II,C_III] = get_fourth_rank_isotropic_tensors
    % ПРИМЕР ИСПОЛЬЗОВАНИЯ:
    % A = randi(9,3,3)
    % 
    % [C_I,C_II] = get_fourth_rank_isotropic_tensors;
    % 
    % C = tensorprod(C_I,A,[3 4],[1 2])
    % % D = tensorprod(C_I,A,4,1)
    % E = tensorprod(C_II,A,[3 4],[1 2])
    % F = tensorprod(A,C_II,[1 2],[1 2])
    % G = tensorprod(C_III,A,[3 4],[1 2])
    % K = tensorprod(A,C_III,[1 2],[1 2])

I = eye(3);

C_I = tensorprod(I,I);

C_II = zeros(3,3,3,3);
for i = 1:3
  for j = 1:3
    for k = 1:3
      for l = 1:3
        C_II(i,j,k,l) = I(i,l) * I(k,j);
      end
    end
  end
end

C_III = zeros(3,3,3,3);
for i = 1:3
  for j = 1:3
    for k = 1:3
      for l = 1:3
        C_III(i,j,k,l) = I(i,k) * I(l,j);
      end
    end
  end
end

end

function [deformation_increments] = get_deformation_increments(...
  time,deformation,fit_model_orientation_case,fit_model_shear_case,...
  sigma,chi_orientation,G,tau_orientation,tau_shear,time_end,...
  theta,alpha,gamma_o,bundles_concentration,kinetic_coefficient_chi)

deformation_increments = zeros(size(deformation));

chi_shear = deformation(3);

dimensionless_orientation_deformation = deformation(1);
dF_dorientation_deformation = theta / gamma_o * (...
  feval(fit_model_orientation_case,dimensionless_orientation_deformation)...
  - chi_orientation^-1 * dimensionless_orientation_deformation... % (lamda_o * gamma_o) / theta * 
  );

dimensionless_shear_deformation = sqrt(alpha / theta) / bundles_concentration * deformation(2);
dF_dshear_deformation = sqrt(alpha * theta) * (...
  feval(fit_model_shear_case,dimensionless_shear_deformation)...
  + dimensionless_shear_deformation...
  - chi_shear.^-1 * dimensionless_shear_deformation...
  );

sigma = get_sigma(sigma,time,time_end);
dimensionless_sigma = sigma / sqrt(alpha * theta); % gamma_o / theta * 

orientation_deformation = -.99:.01:.99;
effective_field_integral_orientation_case = ...
  get_effective_field_integral_orientation_case(orientation_deformation);
F_orientation = ...
  effective_field_integral_orientation_case ...
  - dimensionless_sigma .* orientation_deformation ...
  - chi_orientation^-1 * orientation_deformation.^2 / 2;

F_orientation_barrier = get_energy_barrier(F_orientation * theta);
% F_orientation_barrier = get_current_energy_barrier(F_orientation,orientation_deformation,deformation(1));
tau_orientation = tau_orientation * exp(F_orientation_barrier / (theta));

shear_deformation = 0:.01:5;
effective_field_derivative_shear_case = ...
  get_effective_field_integral_shear_case(shear_deformation);
F_shear = ...
  effective_field_derivative_shear_case  + shear_deformation.^2 / 2 ...
  - dimensionless_sigma .* shear_deformation - chi_shear^-1 * shear_deformation.^2 / 2;

F_shear_barrier = get_energy_barrier(F_shear * theta);
% F_shear_barrier = get_current_energy_barrier(F_shear,shear_deformation,deformation(2));
tau_shear = tau_shear * exp((F_shear_barrier) / theta); % F_orientation_barrier + ?

% if abs(F_shear_barrier) < 1e-1
%   tau = Inf;
% else
%   tau = 1 * exp(F_orientation_barrier / F_shear_barrier);
% end
tau = Inf;

deformation_increments(1) = ...
  1 ./ (G * tau_orientation) .* (sigma - dF_dorientation_deformation) + ...
  1 ./ (G * (tau)) .* (sigma - dF_dshear_deformation);
deformation_increments(2) = 1 ./ (G * (tau_shear)) .* (sigma - dF_dshear_deformation) + ...
  1 ./ (G * (tau)) .* (sigma - dF_dorientation_deformation);
deformation_increments(3) = - kinetic_coefficient_chi * heaviside(deformation(1) - 1) * (2*chi_shear.^2).^-1 .* deformation(2).^2; % 

end

function [sigma] = get_sigma(initial_sigma,time,time_end)

% if time < time_end / 2
  sigma = initial_sigma; %  * time
% else
%   sigma = 0; % initial_sigma * (time_end / 2 - (time - time_end / 2))
% end

end

function [energy_barrier] = get_current_energy_barrier(free_energy,deformation,current_deformation)

[~,n] = min(abs(deformation - current_deformation));

free_energy = free_energy(n:end);

local_maxima = islocalmax(free_energy);
local_minima = islocalmin(free_energy);

if sum(local_maxima) == 0 && sum(local_minima) == 0
  [~,n_maximum] = max(free_energy);
  [~,n_minimum] = min(free_energy);

  if n_maximum == 1
    n_maximum = [];
  elseif n_minimum == 1
    n_minimum = [];
  end

  energy_barrier = ...
    free_energy(min([n_maximum; n_minimum])) - ...
    free_energy(1);
else
  n_first_local_maximum = find(local_maxima,1);
  n_first_local_minimum = find(local_minima,1);
  
  energy_barrier = ...
    free_energy(min([n_first_local_minimum; n_first_local_maximum])) - ...
    free_energy(1);
end

end

function [orientation_deformation_increment] = get_orientation_deformation_increment(...
  time,orientation_deformation,fit_model_orientation_case,sigma,chi_orientation,G,tau)

orientation_deformation_increment = zeros(size(orientation_deformation));

dF_dorientation_deformation = feval(fit_model_orientation_case,orientation_deformation) - orientation_deformation;

deformation = -.99:.01:.99;
effective_field_derivative_orientation_case = ...
  get_effective_field_integral_orientation_case(deformation);
F_orientation = ...
  effective_field_derivative_orientation_case ...
  - sigma * time .* deformation ...
  - chi_orientation^-1 * deformation.^2;

F_orientation_barrier = get_energy_barrier(F_orientation);
tau = tau * exp(F_orientation_barrier);

orientation_deformation_increment(1) = 1 ./ (G * tau) .* (sigma * time - dF_dorientation_deformation);

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

function [effective_field_derivative] = get_effective_field_integral_shear_case(deformation)

effective_field_derivative = ...
  2.97 * atan(.15 * (-5.75 + 2 * deformation)) ...
  - 31.68 * atan(.79 * (2.39 + 2 * deformation)) ...
  + .34 * log(18.82 - 5.75 * deformation + deformation.^2) ...
  + 8.04 * log(1.83 + 2.39 * deformation + deformation.^2) ...
  ;

effective_field_derivative = real(effective_field_derivative);

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

function [effective_field_approximation,fit_model] = get_effective_field_approximation(effective_field,deformation,fit_type,starting_points)

if nargin < 4 || isempty(starting_points)
  fit_model = fit(deformation, effective_field, fit_type);
else
  fit_model = fit(deformation, effective_field, fit_type,'StartPoint',[starting_points(:,1)' starting_points(:,2)']);
end

effective_field_approximation = feval(fit_model,deformation);

end

function [shear_deformation] = get_shear_deformation(effective_field,theta,is_dimensionless,alpha)
% Функция для расчета сдвиговой части деформации, обусловленной
% проскальзыванием пачек актиновых филаментов.

% Проверка входных данных
if nargin < 4 || isempty(alpha)
  alpha = 1;
end

if nargin < 3 || isempty(is_dimensionless)
  is_dimensionless = 0;
end

if nargin < 2 || isempty(theta)
  theta = 1;
end

if nargin < 1 || isempty(effective_field)
  effective_field = (-10:.01:10)';
end

% Расчет сдвиговой части деформации
if is_dimensionless
  shear_deformation = 4 * (...
    - effective_field...
    + ( sqrt(2*pi) * erfi(effective_field / sqrt(8)) ) ./ ...
    hypergeom([.5 .5], [1.5 1.5], effective_field.^2 / 8) ) ./ ...
    effective_field.^2 ...
    ;
else
  shear_deformation = 4 * theta * (...
    - effective_field...
    + ( sqrt(2*pi*alpha*theta) * erfi(effective_field / sqrt(8*alpha*theta)) ) ./ ...
    hypergeom([.5 .5], [1.5 1.5], effective_field.^2 / (8*alpha*theta)) ) ./ ...
    effective_field.^2 ...
    ;
end

end

function [orientation_deformation] = get_orientation_deformation(effective_field,theta,is_dimensionless)
% Функция для расчета ориентационной части деформации, обусловленной
% ориентированием  сегментов актиновых филаментов в направлении сдвига.

% Проверка входных данных
if nargin < 3 || isempty(is_dimensionless)
  is_dimensionless = 0;
end

if nargin < 2 || isempty(theta)
  theta = 1;
end

if nargin < 1 || isempty(effective_field)
  effective_field = (-10:.01:10)';
end

% Расчет ориентационной части деформации
if is_dimensionless
  orientation_deformation = ...
    exp(effective_field / 4) .* sqrt(3/pi * effective_field) ...
    .* (exp(effective_field / 2) ./ effective_field ...
    - (exp(-effective_field / 4) * sqrt(pi/3) ...
    .* (2 + effective_field) .* erfi(1/2 * sqrt(3*effective_field)) ...
    ./ (2*effective_field.^(3/2)))) ./ erfi(1/2 * sqrt(3*effective_field));
else
  orientation_deformation = ...
    exp(effective_field / (4*theta)) .* sqrt(3/pi * effective_field) ...
    .* (exp(effective_field / (2*theta)) * theta ./ effective_field ...
    - (exp(-effective_field / (4*theta)) * sqrt(pi/3 * theta) ...
    .* (2*theta + effective_field) .* erfi(1/2 * sqrt(3*effective_field / theta)) ...
    ./ (2*effective_field.^(3/2)))) / (sqrt(theta) ...
    .* erfi(1/2 * sqrt(3*effective_field / theta)));
end

end

