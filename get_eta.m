function [eta_curve, eta_line] = get_eta(xis,sigma,T,lambda,gamma,R_fun,R_derivative_fun)
% Функция для расчета зависимостей eta(xi), кривой и линии.
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

