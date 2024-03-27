%% Решение уравнения самосогласования
xi_step = .01;
xi = xi_step:xi_step:10;
w =...
  1/6 + 1/4 * sqrt(pi) * xi .* ...
  ( (4 * exp(xi.^2 / 4) + (-xi.^2 ./ 4^(3/4-1)).*igamma(1-3/4,-xi.^2 / 4) ) ./ (4 * xi.^2) + ...
  (sqrt(2) * gamma(5/4)) ./ (-xi.^2).^(5/4) ) + ...
  1/28 * xi.^2 .* hypergeom([1 7/4],[3/2 11/4],xi.^2 / 4);
z = 1/24 * ...
  ( (3*sqrt(2*pi) * (gamma(1/4) - gammainc(-xi.^2/4,1/4))) ./ (-xi.^2).^(1/4)...
  + 4 * xi .* hypergeom([3/4 1],[3/2 7/4],xi.^2 / 4) );
eps_r = real(w ./ z);

xi = xi(eps_r > 0) - min(xi(eps_r > 0));
eps_r = eps_r(eps_r > 0);

%% Аппроксимация решения уравнения самосогласования
p = [5084 3436 4.819e+05 -2216];
q = [178.9 6096 2.403e+05 6.721e+05];
% p = [146.6 6207 3648];
% q = [-65.25  2049 1724];
fit_type = fittype('rat34');
fit_model = fit(eps_r', (xi)', fit_type,'StartPoint',[p q]) % 
fit_xi = feval(fit_model,eps_r);

plot(eps_r,xi)
hold on
plot(eps_r,fit_xi)
% loglog(eps_r,xi)
% hold on
% loglog(eps_r,fit_xi)

%% Влияние структурного параметра
eps_r = 0:.1:30;
chi_c = 3.65;
chi_t = 4.9;
chis = [3 3.6 4 4.5 5];
for chi = chis
  sigma = -((chi.^-1 - 0) .* eps_r - feval(fit_model,eps_r)');

  hold on;
  plot(sigma,eps_r);
end


%% Виды свободной энергии
eps_r = 0:.1:30;
chis = [3.58 4 5];
i = 1;
F_r = zeros(numel(eps_r),3);
for chi = chis
  sigma = .5;
  F_r(:,i) = eps_r.^2 / (2*chi) - (...
    5084 * (-.24 * atan(.014 * (27.35 + 2 * eps_r)) ...
    - .0014 * log(3 + eps_r) + 1.15 * log(148.55 + eps_r) ...
    - .076 * log(1505 + 27.35 * eps_r + eps_r.^2))...
    ) + sigma .* eps_r;
  
  hold on;
  plot(eps_r,-F_r(:,i))
  i = i + 1;
end

%% Моделирование кинетики разрывов актиновых филаментов
chi_c = 3.6;

chi_0 = 3.7;
eps_r_0 = 0;
initial_conditions = [eps_r_0; chi_0];

n_time_steps = 2^12;
time = linspace(0,1.6,n_time_steps);

eps = .1;

C_eps_r = 1;
C_chi = 1;

[time,solution] = ode15s(@(t,previous_step_solution) get_rhs_exprs(t,previous_step_solution,eps,C_eps_r,C_chi),time,initial_conditions);
eps_r = solution(:,1);
chi = solution(:,2);

[rhs_exprs] = get_rhs_exprs(time',solution',eps,C_eps_r,C_chi);
dF_deps_r = rhs_exprs(1,:)';
sigma = 1 * (eps * time - eps_r);
energy_dissipation_rate = eps_r .* (dF_deps_r);

plot(eps*time,energy_dissipation_rate)
hold on
plot(eps*time,(dF_deps_r - C_eps_r * sigma))

pHat = lognfit(diff(eps_r(chi > chi_c)));
sigma_pdf_ductile = pdf('Lognormal',diff(eps_r(chi > chi_c)),pHat(1),pHat(2));
pHat = lognfit(diff(eps_r(chi < chi_c)));
sigma_pdf_brittle = pdf('Lognormal',diff(eps_r(chi < chi_c)),pHat(1),pHat(2));
pHat = lognfit(diff(eps_r));
sigma_pdf = pdf('Lognormal',diff(eps_r),pHat(1),pHat(2));

figure;hold on
plot(time(chi > chi_c)*eps,eps_r(chi > chi_c),'g');
plot(time(chi < chi_c)*eps,eps_r(chi < chi_c),'r');
figure;hold on;
plot(time(chi > chi_c)*eps,chi(chi > chi_c),'g');
plot(time(chi < chi_c)*eps,chi(chi < chi_c),'r');

figure;
plot(diff(eps_r(chi > chi_c)),sigma_pdf_ductile / max(sigma_pdf_ductile),'g'); % 
% loglog(diff(eps_r(chi > chi_c)),sigma_pdf_ductile,'g'); %  / max(sigma_pdf_ductile)
hold on;
plot(diff(eps_r(chi < chi_c)),sigma_pdf_brittle / max(sigma_pdf_brittle),'r'); % 
% loglog(diff(eps_r(chi < chi_c)),sigma_pdf_brittle,'r'); %  / max(sigma_pdf_brittle)
plot(diff(eps_r),sigma_pdf / max(sigma_pdf),'k');
% loglog(diff(eps_r),sigma_pdf / max(sigma_pdf),'k');

%% Моделирование кинетики разрывов актиновых филаментов (НЕПРАВИЛЬНО!)
chi_c = 3.6;

chi_0 = 3.7;
eps_r_0 = 0;
initial_conditions = [eps_r_0; chi_0];

n_time_steps = 2^12;
time = linspace(0,1.6,n_time_steps);

eps = .1;

C_eps_r = 1;
C_chi = 1;

[time,solution] = ode15s(@(t,previous_step_solution) get_rhs_exprs(t,previous_step_solution,eps,C_eps_r,C_chi),time,initial_conditions);
eps_r = solution(:,1);
chi = solution(:,2);

pHat = lognfit(diff(eps_r(chi > chi_c)));
sigma_pdf_ductile = pdf('Lognormal',diff(eps_r(chi > chi_c)),pHat(1),pHat(2));
pHat = lognfit(diff(eps_r(chi < chi_c)));
sigma_pdf_brittle = pdf('Lognormal',diff(eps_r(chi < chi_c)),pHat(1),pHat(2));
pHat = lognfit(diff(eps_r));
sigma_pdf = pdf('Lognormal',diff(eps_r),pHat(1),pHat(2));

figure;hold on
plot(time(chi > chi_c)*eps,eps_r(chi > chi_c),'g');
plot(time(chi < chi_c)*eps,eps_r(chi < chi_c),'r');
figure;hold on;
plot(time(chi > chi_c)*eps,chi(chi > chi_c),'g');
plot(time(chi < chi_c)*eps,chi(chi < chi_c),'r');

figure;
plot(diff(eps_r(chi > chi_c)),sigma_pdf_ductile / max(sigma_pdf_ductile),'g'); % 
% loglog(diff(eps_r(chi > chi_c)),sigma_pdf_ductile,'g'); %  / max(sigma_pdf_ductile)
hold on;
plot(diff(eps_r(chi < chi_c)),sigma_pdf_brittle / max(sigma_pdf_brittle),'r'); % 
% loglog(diff(eps_r(chi < chi_c)),sigma_pdf_brittle,'r'); %  / max(sigma_pdf_brittle)
plot(diff(eps_r),sigma_pdf / max(sigma_pdf),'k');
% loglog(diff(eps_r),sigma_pdf / max(sigma_pdf),'k');


function [solution] = get_rhs_exprs(t,previous_step_solution,eps,C_eps_r,C_chi)

solution = zeros(size(previous_step_solution));
eps_r = previous_step_solution(1,:);
chi = previous_step_solution(2,:);

sigma = 1 * (eps * t - eps_r);

deps_r_dt = C_eps_r * ( (chi.^-1 - 0) .* eps_r - (5084 * ...
  -.0015 ./ (3 + eps_r) ...
  +1.15 ./ (148.5 + eps_r) ...
  -.08 * (27.35 + 2 * eps_r) ./ (1505 + 27.35 * eps_r + eps_r.^2) ...
  -.0067 ./ (1 + .0002 * (27.35 * eps_r + eps_r.^2)) + sigma) ...
   );
dchi_dt = -C_chi ./ (2*chi.^2) .* eps_r.^2;

solution = [deps_r_dt; dchi_dt];

end


