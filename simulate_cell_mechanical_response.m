    time_step = .013;
    time = time_step:time_step:10;

    eta_0 = .04;
    d_eta_0 = .1;
    initial_conditions = eta_0;

    G = 1000;
    dzeta = .0001;
    L_eta = .01;
    L_tau = 0;

%     eps_start = eta_0 * dzeta;
%     eps_end = .1;
%     eps = (eps_start:((eps_end - eps_start)/(length(time)-1)):eps_end);
    eps = .02;
    % deps_dt = gradient(eps) / time_step;

    theta = [.2];
    line_style = {'-','--',':'};

    for i = 1:length(theta)
      [~,eta] = ode15s(@(t,eta) get_rhs_expr(t,eta,eps,G,dzeta,L_tau,L_eta,theta(i)),time,initial_conditions);
      sigma = G * (eps - dzeta * eta(:,1)');

      figure(1);hold on
      plot(eps * time,sigma,'k','LineWidth',2,'LineStyle',line_style{i})
      xlabel(['{\it' char(949) '}']);
      ylabel('\sigma');
      set_figure;

      figure(2);hold on
      plot(eps * time,eta(:,1),'k','LineWidth',2,'LineStyle',line_style{i})
      xlabel(['{\it' char(949) '}']);
      ylabel('\eta');
      set_figure;
    end
    for i = 1:2
      figure(i);legend(compose('{\\it\\Theta}=%.2f',theta));
    end


%%
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

function visualize_double_pl_relaxation_of_cell

time_step = .0013;
time = time_step:time_step:10;
[force_curves] = readmatrix('data.txt')';

n_force_curve = 4; % 1 2

eta_0 = .01;
d_eta_0 = 0.226841041632112; % 1.843893633906291 1.519968801891266
    
L_tau = 1.635154151014720e+03; % 1.038734299154578e+05 7.526189514762888e+04
dzeta = 0.048114417127569; % 0.0054627092184820 0.006464748801383
L_eta = 0.484795734305463; % 0.488093794205235 0.422867214093489
G = 1;
n = 0.014457714548634; % 0.005258019008872 0.004447238153784

eps_start = eta_0 * dzeta;
eps_end = .1;
eps = (eps_start:((eps_end - eps_start)/(length(time)-1)):eps_end);
eps((floor(1*end/10)+1):end) = eps(floor(1*end/10));
    
[~,sigma] = get_force_curve(time,eps,L_tau,dzeta,G,eta_0,d_eta_0,L_eta);

span_alpha = time > .01 & time < 1.;
span_beta = time > 0.0013 & time < .01;
span = time > .0013 & time < 1.;

% x = log(time(span_alpha));
% y = log(force_curves(n_force_curve,span_alpha));
% y_alpha_fitmodel = fit(x',y','(a-0.75)*x + b');
% x = log(time(span_beta));
% y = log(force_curves(n_force_curve,span_beta));
% y_beta_fitmodel = fit(x',y','(a-0.2)*x + b');
% x = log(time(span));
% y = log(force_curves(n_force_curve,span));
% figure(1);
% plot(x,y,'k-')
% hold on
% plot(x,feval(y_alpha_fitmodel,x),'r--')
% plot(x,feval(y_beta_fitmodel,x),'r-.')

x = log(time(span_alpha));
y = log(sigma(span_alpha)/n);
y_alpha_fitmodel = fit(x',y','(a-0.75)*x + b');
x = log(time(span_beta));
y = log(sigma(span_beta)/n);
y_beta_fitmodel = fit(x',y','(a-0.2)*x + b');
x = log(time(span));
y = log(sigma(span)/n);
figure(2);
hold on
plot(x,y,'k-')
plot(x,feval(y_alpha_fitmodel,x),'r--')
plot(x,feval(y_beta_fitmodel,x),'r-.')

% figure(3);hold on
% span = time > 1 & time < 2.9;
% span_beta = time > 1.9 & time < 2.9;
% fit_model = fit(time(span)',sigma(span)'/n,'a+b*exp(c*x)','Start',[1.1 .1 -1]);
% plot(log(time(span)),log(feval(fit_model,time(span)) - fit_model.a))
% fit_model = fit(log(time(span_beta))',log(feval(fit_model,time(span_beta)) - fit_model.a),'a*x+b');
% plot(log(time(span)),feval(fit_model,log(time(span))),'-.')
% span = time > 1 & time < 10;
% span_beta = time > 9 & time < 10;
% fit_model = fit(time(span)',force_curves(n_force_curve,span)','a+b*exp(c*x)','Start',[1 .1 -1]);
% plot(log(time(span)),log(feval(fit_model,time(span)) - fit_model.a))
% fit_model = fit(log(time(span_beta))',log(feval(fit_model,time(span_beta)) - fit_model.a),'a*x+b');
% plot(log(time(span)),feval(fit_model,log(time(span))),'-.')

end

function [n_force_curve] = get_n_force_curve
n_force_curve = 1;
end

function [L_tau, dzeta, L_eta, n] = verificate_parameters(n_force_curve)
% Функция для верификации параметров модели.
%     Пример использования:
%     n_force_curve = 2;  % должно совпадать со значением выставленным в функции get_n_force_curve
%     [L_tau, dzeta, L_eta, d_eta_0, n] = verificate_parameters(n_force_curve);


switch n_force_curve 
  case 1
    coefficients_initial_values = [...
      5 ... % L_tau
      1 ... % dzeta
      1 ... % L_eta
      1]; % n
%       1.843893633906291]; % d_eta_0    

  case 2
    coefficients_initial_values = [...
       4.603758028778943e+02 ... % L_tau
       7.405424894001364 ... % dzeta
       1.181504051673577 ... % L_eta
       2.348331080151025]; % n
%        1.519968801891266]; % d_eta_0

  case 3
    coefficients_initial_values = [...
       350 ... % L_tau
       8 ... % dzeta
       1 ... % L_eta
       1]; % n
%        1.079149000542707]; % d_eta_0

  case 4
    coefficients_initial_values = [...
       1.635154151014720e+03 ... % L_tau
       0.048114417127569... % dzeta
       0.484795734305463 ... % L_eta
       0.014457714548634]; % n
%        0.226841041632112]; % d_eta_0

  case 5
    coefficients_initial_values = [...
       1.677242853312810e+03 ... % L_tau
       0.042351291191231... % dzeta
       0.641116172895703 ... % L_eta
       0.008440505901412]; % n
%        0.260682200799845]; % d_eta_0

  case 6
    coefficients_initial_values = [...
       1.621713960271511e+03 ... % L_tau
       0.036412335929573... % dzeta
       0.964611057733632 ... % L_eta
       0.004242691063008]; % n
%        0.305084842799878]; % d_eta_0
end

% [L_tau, dzeta, L_eta, d_eta_0, n] = identificate_parameters(coefficients_initial_values);
[L_tau, dzeta, L_eta, n] = identificate_parameters(coefficients_initial_values);

end

% function [L_tau, dzeta, L_eta, d_eta_0, n] = identificate_parameters(coefficients_initial_values)
function [L_tau, dzeta, L_eta, n] = identificate_parameters(coefficients_initial_values)
% Функция для идентификации параметров модели L_tau и dzeta.

% проверка входных данных
if sum(size(coefficients_initial_values)) ~= 5
coefficients_initial_values = [...
   7.526189514762888e+04 ... % L_tau
   0.006464748801383... % dzeta
   0.422867214093489 ... % L_eta
   0.004447238153784];... % n
%    1.519968801891266]; % d_eta_0
end

% оптимизация методом Нелдера-Мида
options = optimset('Display','iter','PlotFcns',@optimplotfval);
fun = @get_difference_sigma_and_F;

[coefficients,fval,exitflag,output] = fminsearch(fun,coefficients_initial_values,options)

% выходные данные
L_tau = coefficients(1);
dzeta = coefficients(2);
L_eta = coefficients(3);
n = coefficients(4);
% eta_0 = coefficients(3);
% d_eta_0 = coefficients(5);
end

function [difference_sigma_and_F] = get_difference_sigma_and_F(coefficients)
% Функция для расчета разницы между силовыми кривыми стат. модели и модели КФ.

% загрузка результатов симуляции КФ-модели
time_step = .0013;
time = time_step:time_step:10;
[force_curves] = readmatrix('data.txt')';

% симуляция на основе мезоскопической модели
eta_0 = .01; % .01
% d_eta_0 = coefficients(5); % .9
    
L_tau = coefficients(1);
dzeta = coefficients(2);
L_eta = coefficients(3);
G = 1;

% eps_start = eta_0 * dzeta;
% eps_end = .1;
% eps = (eps_start:((eps_end - eps_start)/(length(time)-1)):eps_end);
% eps((floor(1*end/10)+1):end) = eps(floor(1*end/10));
eps = .01 / time_step;
    
% [~,sigma] = get_force_curve(time,eps,L_tau,dzeta,G,eta_0,d_eta_0,L_eta);
[~,sigma] = get_force_curve(time,eps,L_tau,dzeta,G,eta_0,L_eta);

% диапазон оптимизации
optimization_span = time > 0 & time < 1.1;
n_force_curve = get_n_force_curve;

% сравнение результатов 2 моделей(при оптимизации д.б. закомментировано)
% figure(n_force_curve);hold on;
% plot(time,force_curves(n_force_curve,:),'k-',...
%   time,sigma/coefficients(4),'k--');
% close gcf;

% расчет целевой функции, по которой происходит оптимизация
difference_sigma_and_F = max(abs(sigma(optimization_span)/coefficients(4) - force_curves(n_force_curve,optimization_span)));
end

% function [eta,sigma,time] = get_force_curve(time,eps,L_tau,dzeta,G,eta_0,d_eta_0,L_eta)
function [eta,sigma,time] = get_force_curve(time,eps,L_tau,dzeta,G,eta_0,L_eta)
% Функция для расчета силовой кривой (моделирование экспериметнта на АСМ).
%     Пример использования:
%     time_step = .0013;
%     time = time_step:time_step:10;
%     
%     eta_0 = .01;
%     d_eta_0 = .9;
%     
%     dzeta = .009; % .01 .009 .007 .005
%     
%     eps_start = eta_0 * dzeta;
%     eps_end = .1;
%     eps = (eps_start:((eps_end - eps_start)/(length(time)-1)):eps_end);
%     eps((floor(1*end/10)+1):end) = eps(floor(1*end/10));
%     
%     L_tau = [0 5000 10000]; % 0 5000 10000
%     line_style = {'-','--',':'};
%     
%     for i = 1:length(L_tau)
%       [eta,sigma] = get_force_curve(time,eps,L_tau(i));
%     
%       figure(1);hold on
%       plot(time,sigma,'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('t');
%       ylabel('\sigma');
%       set_figure;
%       
%       figure(2);hold on
%       plot(time,eta(:,1),'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('t');
%       ylabel('\eta');
%       set_figure;
%       
%       figure(3);hold on
%       plot(log(time(2:floor(1*end/10))),log(sigma(2:floor(1*end/10))),'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('t');
%       ylabel('\sigma');
%       set_figure;
%       
%       figure(4);hold on
%       plot(log(time((floor(1*end/10)+1):end)),log(sigma((floor(1*end/10)+1):end)),'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel('t');
%       ylabel('\sigma');
%       set_figure;
%     end
%     for i = 1:4
%       figure(i);legend(compose('{\\itL_{\\tau}}=%i',L_tau));
%     end

if nargin < 8
  L_eta = 1;
end

if nargin < 7
  eta_0 = .01;
  d_eta_0 = .9;
end

if nargin < 5
  G = 1000;
end

if nargin < 4
  dzeta = .009; % .01 .009 .007 .005
end

if nargin < 3
  L_tau = [10000]; % 0 1000 10000
end

if nargin < 2
  eps_start = eta_0 * dzeta;
  eps_end = .1;
  eps = (eps_start:((eps_end - eps_start)/(length(time)-1)):eps_end);
  eps((floor(1*end/10)+1):end) = eps(floor(1*end/10));
end

if nargin < 1
  time_step = .0013;
  time = time_step:time_step:10;
end

% initial_conditions = [eta_0; d_eta_0];
initial_conditions = eta_0;

deps_dt = gradient(eps) / (time(2) - time(1));

% [~,eta] = ode15s(@(t,eta) get_rhs_expr(t,eta,deps_dt,G,dzeta,L_tau,L_eta),time,initial_conditions);
[~,eta] = ode15s(@(t,eta) get_rhs_expr(t,eta,eps,G,dzeta,L_tau,L_eta),time,initial_conditions);
sigma = 3/2 * G * (eps - dzeta * eta(:,1)');

end

function [rhs_expr] = get_rhs_expr(t,varibles,eps,G,dzeta,L_tau,L_eta,theta)
% Функция для расчета правой части дифференциального уравнения
%     Пример использования:
%     time_step = .013;
%     time = time_step:time_step:10;
% 
%     eta_0 = .01;
%     d_eta_0 = .1;
%     initial_conditions = eta_0;
% 
%     G = 1000;
%     dzeta = 1;
%     L_eta = .01;
%     L_tau = 10;
% 
%     eps_start = eta_0 * dzeta;
%     eps_end = .1;
%     eps = (eps_start:((eps_end - eps_start)/(length(time)-1)):eps_end);
%     % deps_dt = gradient(eps) / time_step;
% 
%     theta = [.2 .23 .3];
%     line_style = {'-','--',':'};
% 
%     for i = 1:length(theta)
%       [~,eta] = ode15s(@(t,eta) get_rhs_expr(t,eta,eps,G,dzeta,L_tau,L_eta,theta(i)),time,initial_conditions);
%       sigma = G * (eps - dzeta * eta(:,1)');
% 
%       figure(1);hold on
%       plot(eps,sigma,'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel(['{\it' char(949) '}']);
%       ylabel('\sigma');
%       set_figure;
% 
%       figure(2);hold on
%       plot(eps,eta(:,1),'k','LineWidth',2,'LineStyle',line_style{i})
%       xlabel(['{\it' char(949) '}']);
%       ylabel('\eta');
%       set_figure;
%     end
%     for i = 1:2
%       figure(i);legend(compose('{\\it\\Theta}=%.2f',theta));
%     end


if nargin < 8
  theta = .2;
end

if nargin < 7
  L_eta = 1.;
end

if nargin < 6
  L_tau = 3.0;
end

rhs_expr = zeros(size(varibles));

% rhs_expr(1) = varibles(2);
% rhs_expr(2) = L_eta * (...
%   dzeta * ( G * (eps(ceil(t)) - dzeta * varibles(2)) - L_tau * dzeta * varibles(2) )...
%   + varibles(2)...
%   - theta * varibles(2) * (1.4 / (1 - varibles(1)).^2)...
%   + 2/3 * G * (eps(ceil(t)) - dzeta * varibles(2))...
%   );

sigma = 3/2 * G * (eps * t - dzeta * varibles(1));

eta = -.5:.013:1;
free_energy = get_free_energy_approximation(eta,sigma,theta);
figure(3);hold on;
plot(eta,free_energy);

% dF_deta =...
%   -1.25 / (1 - varibles(1))...
%   -9 / (2 - varibles(1))...
%   -.5 / (.5 + varibles(1))...
%   -1 / (1.1 + varibles(1))...
%   -73.73 / (19.1 + varibles(1))...
%   ;
p1 =        4109;
p2 =       -9339;
p3 =        5731;
p4 =        1115;
p5 =       -1868;
p6 =       132.7;
p7 =       198.2;
p8 =      -26.04;
p9 =       2.659;
p10 =      0.2448;
dF_deta = ...
  p1*varibles(1)^10 ...
  + p2*varibles(1)^9 ...
  + p3*varibles(1)^8 ...
  + p4*varibles(1)^7 ...
  + p5*varibles(1)^6 ...
  + p6*varibles(1)^5 ...
  + p7*varibles(1)^4 ...
  + p8*varibles(1)^3 ...
  + p9*varibles(1)^2 ...
  + p10 * varibles(1);
rhs_expr(1) =...
  (...
  L_eta * dzeta^-1 * ( 2/3 * sigma - 4/9 * theta * dF_deta + varibles(1) )...
  ) /...
  (1 + 2/3 * L_eta * L_tau);

% dF_deta =...
%   -1.25 * varibles(2) / (1 - varibles(1)).^2 ...
%   -9 * varibles(2) / (2 - varibles(1)).^2 ...
%   -.5 * varibles(2) / (.5 - varibles(1)).^2 ...
%   -1 * varibles(2) / (1.1 - varibles(1)).^2 ...
%   -73.73 * varibles(2) / (19.1 - varibles(1)).^2 ...
%   ;
% rhs_expr(1) = varibles(2);
% rhs_expr(2) = L_eta * (...
%   G * eps(ceil(t))...
%   + varibles(2)...
%   - 4/9 * theta * dF_deta ... % varibles(2) * (1.4 / (1 - varibles(1)).^2)
%   ) / ...
%   (dzeta * (1 + 2/3 * L_eta * L_tau));

end

function [dF_deta_uniaxial_ordering] = get_dF_deta_uniaxial_ordering(eta, theta, Sigma)
% Функция для расчета частной производной свободной энергии по eta
% (случай аппроксимации натурального логарифма интеграла в виде a*exp(b*x) + c*exp(d*x)).
%     Пример использования:
%     theta = .2;
%     Sigma = -.15; % .2 .15 -.15
%     eta = -.5:.01:1;
%     
%     [dF_deta_uniaxial_ordering] = get_dF_deta_uniaxial_ordering(eta, theta, Sigma);
%     
%     figure(1);hold on;
%     plot(eta,dF_deta_uniaxial_ordering,'k');
%     xlabel('\eta');
%     ylabel('dF/d\eta');
% %     hold off;


etas = -.5:.01:1;
xis = 9 / (4 * theta) * (etas + 2/3 * Sigma);
[~,~,fitmodel] = get_R_approximations(xis);

a = fitmodel.a;
b = fitmodel.b;
c = fitmodel.c;
d = fitmodel.d;

dF_deta_uniaxial_ordering = 3/2 * eta + 3/4 - theta...
  * (...
    (9 * a * (b - d) * exp(9 * b * (eta+(2 * Sigma) / 3) / (4 * theta)))...
    ./ (4 * theta * (a * exp(9 * b * (eta+(2 * Sigma) / 3) / (4 * theta)) ...
        + c * exp(9 * d * (eta+(2 * Sigma) / 3) / (4 * theta))))...
    + (9 * d) / (4 * theta)...
    );

end        

function [fitobj] = get_xi_approximation(xis)
% Функция для расчета вида аппроксимации обратной функции xi.
% General model:
%      fitobj(x) = (  a1 * x^4 + b1 * x^3 + c1 * x^2 + d1 * x^1 + e1 ) / ( a2 * 
%                     x^5 + b2 * x^4 + c2 * x^3 + d2 * x^2 + e2 * x^1 + f2 )
%      Coefficients (with 95% confidence bounds):
%        a1 =     -0.0115  (-1.589, 1.566)
%        a2 =   0.0001344  (-0.1453, 0.1456)
%        b1 =     -0.0169  (-3.171, 3.137)
%        b2 =    0.002377  (-0.2986, 0.3034)
%        c1 =     0.02019  (-1.622, 1.662)
%        c2 =   -0.003805  (-0.4592, 0.4516)
%        d1 =     0.01832  (-2.446, 2.483)
%        d2 =   -0.005388  (-0.9529, 0.9421)
%        e1 =  -0.0001432  (-0.01897, 0.01868)
%        e2 =    0.003987  (-0.3399, 0.3479)
%        f2 =    0.002682  (-0.3559, 0.3613)

if nargin < 1
  xis = -20:.13:20; % -4:.13:12
end

syms xi x;

R_fun = exp(xi .* x.^2);
R_fun_der = (x.^2) .* exp(xi .* x.^2);

int_R_fun = int(R_fun,x,[0 1]);
int_R_fun = (pi^(1/2)*erfi((xi)^(1/2)))/(2*(xi)^(1/2));
int_R_fun_der = int(R_fun_der,x,[0 1]);

eta_fun = matlabFunction(3/2 * (int_R_fun_der/int_R_fun - 1/3));

eta = eta_fun(xis);

fitobj = fit(eta', xis',...
  '(  a1 * x^4 + b1 * x^3 + c1 * x^2 + d1 * x^1 + e1 ) / ( a2 * x^5 + b2 * x^4 + c2 * x^3 + d2 * x^2 + e2 * x^1 + f2 )');
% fitobj = fit(eta', xis',...
%   '(  d1 * x^1 + e1 ) / ( d2 * x^2 + e2 * x^1 + f2 )')
end

function [R_approximation,R,fitmodel] = get_R_approximations(xi)
% Функция для расчета аппроксимации R(xi).
%     Пример использования:
%     xi = -5:.06:5;
% 
%     [R_approximation,R] = get_R_approximations(xi);
% 
%     figure(1);hold on;
%     plot(xi, R,'k');
%     plot(xi, R_approximation,'r');
%     xlabel('\xi')
%     ylabel('R');
%     hold off;


R = (sqrt(pi) * erfi(sqrt(xi))) ./ (2 * sqrt(xi));
R = R(~isnan(R));
xi = xi(~isnan(R));
R = R(~isinf(R));
xi = xi(~isinf(R));


fitmodel = fit(xi', R','exp2','Lower',[0,0],'Upper',[3,3],'StartPoint',[.1 .1 .1 .1]); % a*exp(b*x) + c*exp(d*x))
% fitmodel = fit(xi', R', '1 / ( d1 * x^2 + e1 * x^1 + f1 ) + 1 / ( d2 * x^2 + e2 * x^1 + f2 )');
R_approximation = feval(fitmodel,xi);

end

function [R_form] = get_R_form(xi,x)
% Функция для получения вида R(xi).
% В Wolfram Mathematica R(xi) = (sqrt(pi) * erfi(sqrt(xi))) ./ (2 * sqrt(xi))

if nargin < 2
  syms xi x;
end

R_fun = exp(xi .* x.^2);
R_form = matlabFunction(int(R_fun,x,[0 1]));

end

function [free_energy_derivative] = get_free_energy_approximation_derivative(eta,sigma,T,lambda,gamma)
% Функция для расчета производной по eta аппроксимации свободной энергии.

if nargin < 5
  lambda = 1;
  gamma = 1;
end

if nargin < 3
  T = .2;
end

if nargin < 2
  sigma = .2; % -.15 .15 .2
end

if nargin < 1
  eta = -.5:.01:1;
end

% a1 =    1.824e+08; % -0.0115
% a2 =    -4.19e+07; % 0.0001344
% b1 =   -1.252e+08; % -0.0169
% b2 =     6.07e+07; % 0.002377
% c1 =   -1.115e+08; % 0.02019
% c2 =    3.019e+07; % -0.003805
% d1 =    7.913e+06; % 0.01832
% d2 =   -3.875e+07; % -0.005388
% e1 =    1.041e+06; % -0.0001432
% e2 =   -1.237e+07; % 0.003987
% f2 =     2.17e+06; % 0.002682
% P = ( a1 * eta.^4 + b1 * eta.^3 + c1 * eta.^2 + d1 * eta.^1 + e1 );
% M = ( a2 * eta.^5 + b2 * eta.^4 + c2 * eta.^3 + d2 * eta.^2 + e2 * eta.^1 + f2 );
% xi_approximation_derivative = P ./ M;

% xi_approximation_derivative =...
%   1.25 ./ (1 - eta)...
%   + 9 ./ (2 - eta)...
%   - .5 ./ (.5 + eta)...
%   - 1. ./ (1.1 + eta)...
%   - 73.73 ./ (19.1 + eta)...
%   ;

% xi_approximation_derivative =...
%   (2.576e+06 * eta + 7.964e+04) ./ ( -1.34e+06 * eta.^2 + 6.695e+05 * eta + 6.661e+05 );

xi_approximation_derivative = 1.40034 ./ (0.997369 - 1. * eta);

free_energy_derivative = ...
  - eta ...
  + (4 * T) / (9 * lambda) * xi_approximation_derivative...
  - (2 * gamma * sigma) / (3 * lambda)...
  ;

end

function [free_energy_uniaxial_ordering] = get_free_energy_approximation(eta,sigma,T,lambda,gamma)
% Функция для расчета аппроксимации свободной энергии.
%     Пример использования:
%     eta = -.5:.013:1;
%     free_energy = get_free_energy_approximation(eta);
%     
%     figure(1);hold on;
%     plot(eta,free_energy);
%     xlabel('\eta');
%     ylabel('F^o');
%     hold off;

if nargin < 5
  lambda = 1;
  gamma = 1;
end

if nargin < 3
  T = .2;
end

if nargin < 2
  sigma = .2; % -.15 .15 .2
end

if nargin < 1
  eta = -.5:.013:1;
end

xi_approximation =...
  - 1.25 * log(1 - eta)...
  - 9 * log(2 - eta)...
  - .5 * log(.5 + eta)...
  - log(1.1 + eta)...
  - 73.73 * log(19.1 + eta);
% xi_approximation =...
%   -1.40034 * log(0.997369 - 1. * eta); % - 0.267037 * log(0.498328 + 1. * eta)

free_energy_uniaxial_ordering = ...
  - .5 * eta.^2 ...
  + (4 * T) / (9 * lambda) * xi_approximation...
  - (2 * gamma * sigma) / (3 * lambda) * eta...
  + 19.8085; % + 19.8085

end

function [free_energy_uniaxial_ordering] = get_free_energy_uniaxial_ordering(eta, theta, Sigma)
% Функция для расчета свободной энергии одноосного упорядочения полимеров.
%     Пример использования:
%     theta = .2;
%     Sigma = .15; % .2 .15 -.15
%     eta = -.5:.01:1;
% 
%     [free_energy_uniaxial_ordering] = get_free_energy_uniaxial_ordering(eta, theta, Sigma);
% 
%     figure(1);hold on;
%     plot(eta, free_energy_uniaxial_ordering, 'k');
%     xlabel('\eta');
%     ylabel('F/n\lambda');
%     % hold off;

xis = 9 / (4 * theta) * (eta + 2/3 * Sigma);

% Решение с интегралом
R_fun = @(x,xi) exp(xi .* x.^2);
R = zeros(size(eta));
i = 1;
for xi = xis
  R(i) = integral(@(x) R_fun(x,xi),0,1);
  i = i + 1;
end

free_energy_uniaxial_ordering = 3/4 * eta .* (eta + 1)...
  + 1/2 * Sigma - theta * log(R);

% Решение с взятым интегралом и его аппроксимацией
% [R_app, R] = get_R_approximations(xis);
% free_energy_uniaxial_ordering = 3/4 * eta .* (eta + 1)...
%   + 1/2 * Sigma - theta * log(R); % R R_app'

end

function [Sigma, eta] = get_eta_Sigma(T, sigmas)
% Функция для расчета зависимости eta(Sigma).
%     Пример использования:
%     sigmas = -.1:.01:.2;
%     T = 0.3; % .3 .32 .346 .4
% 
%     [Sigma, eta] = get_eta_Sigma(T, sigmas);
% 
%     figure(1);hold on;
%     plot(Sigma, eta, 'k.');
%     xlabel('\Pi');
%     ylabel('\eta');
%     hold off;

if nargin < 2
  sigmas = -.2:.01:.2;
end

if nargin < 1
  T = 0.03;
end

xis = -10:.1:10;

i = 1;
eta = [];
Sigma = [];
for s = sigmas

  [eta_curve, eta_line] = get_eta(xis,s,T);

  points = find(sign((eta_curve - eta_line')) ~= sign(circshift(eta_curve - eta_line',1)));
  points = points(2:end);
  
  for j = 1:length(points)
    eta(i) = eta_curve(points(j));
    Sigma(i) = s;
    i = i + 1;
  end
end

end

function [theta, eta] = get_eta_theta(sigma, Ts)
% Функция для расчета зависимости eta(T).
%     Пример использования:
%     sigma = .0; % 0 .01 .02 .2 -.02 -.2
%     Ts = 0.01:.01:.4;
% 
%     [theta, eta] = get_eta_theta(sigma, Ts);
% 
%     figure(1);hold on;
%     plot(theta, eta, 'k.');
%     xlabel('\Theta');
%     ylabel('\eta');
%     hold off;

if nargin < 2
  Ts = 0.03:.001:.08;
end

if nargin < 1
  sigma = 0;
end

xis = -30:.5:30;

i = 1;
eta = 0;
theta = 0;
for T = Ts

  [eta_curve, eta_line] = get_eta(xis,sigma,T);

  points = find(sign((eta_curve - eta_line')) ~= sign(circshift(eta_curve - eta_line',1)));
  points = points(2:end);
  
  for j = 1:length(points)
    eta(i) = eta_curve(points(j));
    theta(i) = T;
    i = i + 1;
  end
end

end

function [eta_curve, eta_line] = get_eta(xis,sigma,T,lambda,gamma,R_fun,R_derivative_fun)
% Функция для расчета зависимости eta(xi).
%     Пример использования:
%     T = .4; % .4 .26 .16
%     sigma = .1; % .1 -.5 .5
% 
%     xis = -10:.1:20;
% 
%     [eta_curve, eta_line] = get_eta(xis,sigma,T);
% 
%     figure(1);hold on;
%     plot(xis, eta_curve, 'k');
%     plot(xis, eta_line, 'k');
%     ylim([min(eta_curve) max(eta_curve)]);
%     xlabel('\xi');
%     ylabel('\eta');
%     hold off;
%     set_figure;


if nargin < 7
  R_fun = @(x,xi) exp(xi .* x.^2);
  R_derivative_fun = @(x,xi) (x.^2) .* exp(xi .* x.^2);
end

if nargin < 4
  lambda = 1; % .1 eV   .68
  gamma = 1; % 1e-23 sm^3
end

if nargin < 2
  T = .07;
  sigma = 0;
end

R = zeros(length(xis),1);
R_derivative = zeros(size(R));
i = 1;
for xi = xis
  R(i) = integral(@(x) R_fun(x,xi),0,1);
  R_derivative(i) = integral(@(x) R_derivative_fun(x,xi),0,1);
  i = i + 1;
end

eta_curve = 3/2 * (R_derivative ./ R - 1/3);
eta_line = 4*T / (9*lambda) * xis - (2*gamma*sigma) / (3*lambda);

end

function set_figure
% Функция для настройки внешнего вида графиков.

% Font settings
font = 'Times New Roman';
fontSize = 14;
fontWeight = 'normal';

box on;

% Labels
%{
xLabel = '№ гена';
xlabel(xLabel);
yLabel = 'Q''';
ylabel(yLabel);
%}

% Axis and grid settings
%{
grid off;
axis xy;
axis tight;
%}

% Axis limits
%{
% xlim([0 600]);
ylim([0.1 0.17]);
%}

%{
colormap('jet');
c = colorbar;
c.Label.String = 'log_{10}(|A(u,v)|)';
%}

% Title
% title('MCF-7 with BMK');

set(gca,'FontName',font,'FontSize',fontSize,'FontWeight',fontWeight);
end

