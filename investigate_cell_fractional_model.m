E = 1000;
tau = 1e-1;
alpha = .75;
beta = .2;
lambda = 2;

time_step = .0013;
time = (time_step:time_step:1)';

indentation = linspace(0,1,length(time))';

force = E * gamma(lambda + 1) * (time / time(end)).^lambda ...
  .* ( 1 / gamma(lambda + 1 - alpha) * (time / tau).^-alpha ...
      + 1 / gamma(lambda + 1 - beta) * (time / tau).^-beta );

relaxation_function_frac_model = E / gamma(1 - alpha) * (time / tau).^-alpha ...
  + E / gamma(1 - beta) * (time / tau).^-beta;
[relaxation_young_modulus] = get_relaxation_young_modulus(force,indentation,time);

n_scale = length(time) - 1;
scale_size = 1e-2 * time(end);
max_scale = 2^6;
[wt_coefficients_0_order,wt_coefficients_2_order] = get_wavelet_transform_coefficients(force(1:end-1),time(1:end-1),scale_size,n_scale,max_scale);

[F] = get_sigma_with_relaxation_young_modulus(relaxation_function_frac_model,indentation,time);

figure(1);hold on;
% yyaxis left;
plot((time),(relaxation_function_frac_model));
% plot(log10(time),log10(relaxation_young_modulus));
% plot((time(1:end-1)),(max(wt_coefficients_2_order,[],1)));
% yyaxis right;
% plot(time,force);
% plot(time(1:end-1),wt_coefficients_0_order(:,1))
% plot(time,F)
% xline(0)


