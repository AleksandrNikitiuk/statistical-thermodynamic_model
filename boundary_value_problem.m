function [eta] = get_eta_from(nonlocality_parameters,eps,theta,time)
% Функция для исследования зависимости Y от параметра нелокальности.
%     Пример использования:
%     time = linspace(0,(2^4 - 1),2^10); % 2^9 2^11 2^10
%     dt = time(2) - time(1);
%     eps = 0:(.05/(length(time)-1)):.05;
%     theta = .2;
%     nonlocality_parameters = [0.]; % 0 .1 .2
%     [eta] = get_eta_from(nonlocality_parameters,eps,theta,time);
% 
%     model_parameters = get_model_parameters;
%     % figure; hold on;
%     % for i = 1:length(nonlocality_parameters)
%     %   plot(1:model_parameters.n_elements,eta{i}(end,:));
%     % end
%     % hold off;
%     % legend(compose('{\\itA}''=%.1f',nonlocality_parameters));
%     % set_figure;
%     % xlabel('№ element');
%     % ylabel('{\it\eta}');
%     figure(1);
%     hold on;
%     mesh(model_parameters.G * (eps' - model_parameters.dzeta * eta{1}));
%     % plot(eps,sum(model_parameters.G * (eps' - model_parameters.dzeta * eta{1}),2));
%     % plot(eps,model_parameters.G * (eps' - model_parameters.dzeta * eta{1}(:,64)));
%     % plot(eta{1}(:,100));


eta = cell(1);
symmetry_type = 0; % 0 - декартова система координат
n_nonlocality_parameters = length(nonlocality_parameters);
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
model_parameters = get_model_parameters;

for n_nonlocality_parameter = 1:n_nonlocality_parameters
  A = nonlocality_parameters(n_nonlocality_parameter);
  sol = pdepe(...
    symmetry_type,...
    @(z,time,eta,deta_dz) get_pde_terms_nonlocality(z,time,eta,deta_dz,A,eps,theta),...
    @get_ic_nonlocality,...
    @get_bc_nonlocality,...
    model_parameters.elements,...
    time,...
    options);

  eta{n_nonlocality_parameter,1} = sol(:,:,1);
end

end

function [c,f,s] = get_pde_terms_nonlocality(z,time,eta,deta_dz,A,eps,theta)
% Функция для получения членов и коэффициентов уравнения для eta в частных
% производных с учетом эффектов нелокальности.

model_parameters = get_model_parameters;
nonlocality_term = A * deta_dz;

c = 1 - model_parameters.L_eta * model_parameters.L_tau * model_parameters.dzeta^2;
f = model_parameters.L_eta * (nonlocality_term); 
s = model_parameters.L_eta * (...
  (model_parameters.dzeta - 2/3) * model_parameters.G * (eps(floor(time)+1) - model_parameters.dzeta * eta)...
  - theta ...
  * (1.25 ./ (1 - eta)...
  + 9 ./ (2 - eta)...
  - .5 ./ (.5 + eta)...
  - 1. ./ (1.1 + eta)...
  - 73.73 ./ (19.1 + eta))...
  - eta...
  );
end

function [eta_0] = get_ic_nonlocality( z )
% Функция для получения начальных условий уравнения для eta в частных
% производных с учетом эффектов нелокальности.

if z == 2^6
  eta_0 = .0;
else 
  eta_0 = .0;
end

model_parameters = get_model_parameters;
alpha = .29;
eta_0 = ((-1).^(z) .* 1e-0) ./...
  cosh(alpha * ( z - (model_parameters.n_elements / 2)));

end

function [pl,ql,pr,qr] = get_bc_nonlocality(zl,eta_l,zr,eta_r,time)
% Функция для получения граничных условий уравнения для eta в частных
% производных с учетом эффектов нелокальности, которые формируются
% согласно соотношению p + q * f(time, z, Q, dQ_dz) = 0.

% Периодические граничные условия
% pl = eta_l - eta_r;
% ql = 0;
% pr = eta_r - eta_l;
% qr = 0;

% Граниченые условия 1го рода
% pl = eta_l - .0;
% ql = 0;
% pr = eta_r + .0;
% qr = 0;

% Граниченые условия 2го рода
pl = 0.0;
ql = 1;
pr = 0.0;
qr = 1;

% Граниченые условия 3го рода
% pl = -1e-12 * eta_l;
% ql = 1;
% pr = -1e-12 * eta_r;
% qr = 1;
end

function [model_parameters] = get_model_parameters
% Функция для получения параметров модели.

model_parameters.d_eta_0 = 0.1; % 1.843893633906291 1.519968801891266
model_parameters.L_tau = 1; % 1.038734299154578e+05 7.526189514762888e+04
model_parameters.dzeta = .043; % 0.0054627092184820 0.006464748801383
model_parameters.L_eta = 1; % 0.488093794205235 0.422867214093489
model_parameters.G = 100;
model_parameters.n = 1; % 0.005258019008872 0.004447238153784
model_parameters.n_elements = 2^7;
model_parameters.elements = linspace(1,model_parameters.n_elements,model_parameters.n_elements);
end