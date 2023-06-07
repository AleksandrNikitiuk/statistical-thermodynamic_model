function [sigma, G] = get_sigma_frac_model(eps,deps_dt,time,alpha,beta,E,tau)

for t = 1:length(time)
  sigma(t) = trapz()
end

end

function relaxation_function = get_relaxation_function(time,alpha,beta,E,tau)

relaxation_function = ...
  E * (1/gamma(1 - alpha) * (time/tau)^(-alpha) ...
     + 1/gamma(1 - beta) * (time/tau)^(-beta));

end