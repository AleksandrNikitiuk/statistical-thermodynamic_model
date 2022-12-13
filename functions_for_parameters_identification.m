function [n_force_curve] = get_n_force_curve
n_force_curve = 1;
end

function [E2, Lo, nu, n] = verificate_parameters(n_force_curve)
% Функция для верификации параметров модели.
%     Пример использования:
%     n_force_curve = 1;  % должно совпадать со значением выставленным в функции get_n_force_curve
%     [E2, Lo, nu, n] = verificate_parameters(n_force_curve);


switch n_force_curve 
  case 1
    coefficients_initial_values = [...
      100 ... % E2 30
      1 ... % Lo 1
      100 ... % nu .05
      1 ]; % n .001

  case 2    

  case 3    

  case 4    

  case 5    

  case 6
    
end

[E2, Lo, nu, n] = identificate_parameters(coefficients_initial_values);

end

function [E2, Lo, nu, n] = identificate_parameters(coefficients_initial_values)
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
fun = @get_difference_f_and_F;

[coefficients,fval,exitflag,output] = fminsearch(fun,coefficients_initial_values,options)

% выходные данные
E2 = coefficients(1);
Lo = coefficients(2);
nu = coefficients(3);
n = coefficients(4);
end

function [difference_f_and_F] = get_difference_f_and_F(coefficients)
% Функция для расчета разницы между силовыми кривыми стат. модели и модели КФ.

% загрузка результатов симуляции КФ-модели
time_step = .0013;
time = (time_step:time_step:10)';
[force_curves] = readmatrix('data.txt')';

% симуляция на основе мезоскопической модели
time_loading = time(time <= 1);
time_dwell = time(time > 1);

sigma_0 = .0;
    
E1 = 300;
E2 = coefficients(1);
Lo = coefficients(2);
nu = coefficients(3);

eps = .005;
deps_dt = eps * time_step;

theta = .2;
    
[eps_t,sigma] = get_force_curve(time_loading,time_dwell,sigma_0,eps,deps_dt,E1,E2,Lo,nu,theta);
f = (sigma ./ eps_t) .* eps_t.^(3/2);

% диапазон оптимизации
optimization_span = time > 0 & time < 11;
n_force_curve = get_n_force_curve;

% сравнение результатов 2 моделей(при оптимизации д.б. закомментировано)
figure(n_force_curve);hold on;
plot(time,force_curves(n_force_curve,:),'k-',...
  time,f/coefficients(4),'k--');
close gcf;

% расчет целевой функции, по которой происходит оптимизация
difference_f_and_F = max(abs(f(optimization_span)/coefficients(4) - force_curves(n_force_curve,optimization_span)'));
end

function [eps_t,sigma,time] = get_force_curve(time_loading,time_dwell,sigma_0,eps,deps_dt,E1,E2,Lo,nu,theta)
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

if nargin < 13
  lambda = 1.;
  gamma = 1.;
end

if nargin < 10
  theta = .2;
end

if nargin < 9
  nu = 1.;
end

if nargin < 8
  Lo = 1.;
end

if nargin < 7
  E2 = 1.;
end

if nargin < 6
  E1 = 1.;
end

if nargin < 5
  deps_dt = 0.;
end

if nargin < 4
  eps = .05;
end

if nargin < 3
  sigma_0 = 0;
end

if nargin < 2
  time_step = .013;
  time = (time_step:time_step:10)';
  time_loading = time(time <= 1);
  time_dwell = time(time > 1);
end

initial_conditions = sigma_0;

[~,sigma] = ode15s(@(t,sigma) get_rhs_expr(t,sigma,eps,deps_dt,E1,E2,Lo,nu,theta),time_loading,initial_conditions);
initial_conditions = sigma(end);
[~,sigma_dwell] = ode15s(@(t,sigma) get_rhs_expr(t,sigma,eps,0,E1,E2,Lo,nu,theta),time_dwell,initial_conditions);
sigma = [sigma; sigma_dwell];

eps_t = [eps*time_loading; eps * ones(length(time_dwell),1)];
time = [time_loading; time_dwell];

end