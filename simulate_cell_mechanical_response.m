


%%
function [eta_curve, eta_line] = get_eta(xis,sigma,T,lambda,gamma,R_fun,R_derivative_fun)
% ������� ��� ������� ����������� eta(xi).
%     ������ �������������:
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

function [free_energy_uniaxial_ordering] = get_free_energy_uniaxial_ordering(eta, theta, Sigma)
% ������� ��� ������� ��������� ������� ���������� ������������ ���������.
%     ������ �������������:
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

% ������� � ����������
R_fun = @(x,xi) exp(xi .* x.^2);
R = zeros(size(eta));
i = 1;
for xi = xis
  R(i) = integral(@(x) R_fun(x,xi),0,1);
  i = i + 1;
end

free_energy_uniaxial_ordering = 3/4 * eta .* (eta + 1)...
  + 1/2 * Sigma - theta * log(R);

% ������� � ������ ���������� � ��� ��������������
% [R_app, R] = get_R_approximations(xis);
% free_energy_uniaxial_ordering = 3/4 * eta .* (eta + 1)...
%   + 1/2 * Sigma - theta * log(R); % R R_app'

end

function [dF_deta_uniaxial_ordering] = get_dF_deta_uniaxial_ordering(eta, theta, Sigma)
% ������� ��� ������� ������� ����������� ��������� ������� �� eta
% (������ ������������� ������������ ��������� ��������� � ���� a*exp(b*x) + c*exp(d*x)).
%     ������ �������������:
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

function [R_approximation,R,fitmodel] = get_R_approximations(xi)
% ������� ��� ������� ������������� R(xi).
%     ������ �������������:
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
% ������� ��� ��������� ���� R(xi).
% � Wolfram Mathematica R(xi) = (sqrt(pi) * erfi(sqrt(xi))) ./ (2 * sqrt(xi))

if nargin < 2
  syms xi x;
end

R_fun = exp(xi .* x.^2);
R_form = matlabFunction(int(R_fun,x,[0 1]));

end

function [free_energy_derivative] = get_free_energy_approximation_derivative(eta,sigma,T,lambda,gamma)
% ������� ��� ������� ����������� �� eta ������������� ��������� �������.

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
% ������� ��� ������� ������������� ��������� �������.
%     ������ �������������:
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
% ������� ��� ������� ����������� eta(Sigma).
%     ������ �������������:
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
% ������� ��� ������� ����������� eta(T).
%     ������ �������������:
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

function set_figure
% ������� ��� ��������� �������� ���� ��������.

% Font settings
font = 'Times New Roman';
fontSize = 14;
fontWeight = 'normal';

box on;

% Labels
%{
xLabel = '� ����';
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

function [eta] = get_eta_from(nonlocality_parameters,eps,theta,time)
% ������� ��� ������������ ����������� Y �� ��������� �������������.
%     ������ �������������:
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
%     % xlabel('� element');
%     % ylabel('{\it\eta}');
%     figure(1);
%     hold on;
%     mesh(model_parameters.G * (eps' - model_parameters.dzeta * eta{1}));
%     % plot(eps,sum(model_parameters.G * (eps' - model_parameters.dzeta * eta{1}),2));
%     % plot(eps,model_parameters.G * (eps' - model_parameters.dzeta * eta{1}(:,64)));
%     % plot(eta{1}(:,100));


eta = cell(1);
symmetry_type = 0; % 0 - ��������� ������� ���������
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
% ������� ��� ��������� ������ � ������������� ��������� ��� eta � �������
% ����������� � ������ �������� �������������.

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
% ������� ��� ��������� ��������� ������� ��������� ��� eta � �������
% ����������� � ������ �������� �������������.

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
% ������� ��� ��������� ��������� ������� ��������� ��� eta � �������
% ����������� � ������ �������� �������������, ������� �����������
% �������� ����������� p + q * f(time, z, Q, dQ_dz) = 0.

% ������������� ��������� �������
% pl = eta_l - eta_r;
% ql = 0;
% pr = eta_r - eta_l;
% qr = 0;

% ���������� ������� 1�� ����
% pl = eta_l - .0;
% ql = 0;
% pr = eta_r + .0;
% qr = 0;

% ���������� ������� 2�� ����
pl = 0.0;
ql = 1;
pr = 0.0;
qr = 1;

% ���������� ������� 3�� ����
% pl = -1e-12 * eta_l;
% ql = 1;
% pr = -1e-12 * eta_r;
% qr = 1;
end

function [model_parameters] = get_model_parameters
% ������� ��� ��������� ���������� ������.

model_parameters.d_eta_0 = 0.1; % 1.843893633906291 1.519968801891266
model_parameters.L_tau = 1; % 1.038734299154578e+05 7.526189514762888e+04
model_parameters.dzeta = .043; % 0.0054627092184820 0.006464748801383
model_parameters.L_eta = 1; % 0.488093794205235 0.422867214093489
model_parameters.G = 100;
model_parameters.n = 1; % 0.005258019008872 0.004447238153784
model_parameters.n_elements = 2^7;
model_parameters.elements = linspace(1,model_parameters.n_elements,model_parameters.n_elements);
end

