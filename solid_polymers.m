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

