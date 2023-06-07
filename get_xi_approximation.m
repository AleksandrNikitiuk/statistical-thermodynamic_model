function [fit_xis,fit_model] = get_xi_approximation(xis,type_approximation)
% Функция для расчета вида аппроксимации обратной функции xi.
%     Пример использования:
%     xis = -10:.1:20;
% 
%     [fit_xis,fit_model] = get_xi_approximation(xis);
%     eta_curve = get_eta(xis);
% 
%     figure(1);hold on;
%     plot(eta_curve, xis, 'k');
%     % plot(eta_curve, feval(fit_model,eta_curve))
%     plot(eta_curve, fit_xis)
%     xlabel(['{\it' char(949) '}_o'])
%     ylabel('\xi')
% 
%     fit_model


if nargin < 1
  xis = -10:.1:20;
end

if nargin < 2
  type_approximation = 2; % 1 - сумма эспонент, 2 - отношение полиномов
end

eta_curve = get_eta(xis);

switch type_approximation
  case 1
    a =    -0.01655;%  (-0.01918, -0.01391)
    b =      -14.07;%  (-14.46, -13.68)
    c =   1.189e-05;%  (9.109e-06, 1.468e-05)
    d =       15.04;%  (14.79, 15.29)
    e =       7.506;%  (7.419, 7.593)
    fit_options = fitoptions('Method','NonlinearLeastSquares',...
      'StartPoint',[a b c d e]);
    fit_type = fittype('a*exp(b*x)+c*exp(d*x)+e*x','options',fit_options); % 
    fit_model = fit(eta_curve, xis', fit_type);
    fit_xis = a*exp(b*eta_curve)+c*exp(d*eta_curve)+e*eta_curve;
  case 2
    p1 =  -1.525e+05;
    p2 =   3.533e+05;
    p3 =   1.452e+04;
    q1 =  -8.702e+04;
    q2 =   3.854e+04;
    q3 =   4.916e+04;
    fit_type = fittype('rat23');
    fit_model = fit(eta_curve, xis', fit_type,'StartPoint',[p1, p2, p3, q1, q2, q3]);
    fit_xis = feval(fit_model,eta_curve);
end

end