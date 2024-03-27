function [relaxation_young_modulus] = get_relaxation_young_modulus(force,indentation,time,beta,alpha)
% Функция для расчета релаксационного модуля Юнга по данным АСМ.
%     Пример использования:
%     time_step = .013;
%     time = (time_step:time_step:10)';
%     time_loading = time(time <= 1);
%     
%     Fs = readmatrix('data.txt');
%     F = Fs(:,1);
%     F_loading = F(time <= 1);
%     
%     [relaxation_young_modulus] = get_relaxation_young_modulus(F_loading,[],time_loading,3/2);
%     plot(time_loading, relaxation_young_modulus);



% Проверка входных данных
if nargin < 1 || isempty(force)
  error('Please pass force values to the function.');
end

n = length(force);
if sum(size(force)) ~= (n + 1)
  error('The force must be a vector (1D-array).');
end

if nargin < 5 || isempty(alpha)
  alpha = 1;
end

if nargin < 4 || isempty(beta)
  beta = 3/2;
end

if nargin < 3 || isempty(time)
  time = 0:(n - 1);
end

if nargin < 2 || isempty(indentation)
  indentation = zeros(size(force));
  indentation(:,:) = 1 / (n - 1):(1 - 1 / (n - 1)) / (n - 1):1;
end

% Расчет релаксационного модуля Юнга на основе метода трапеций с const шагом по времени 
time_step = time(2) - time(1);

force_increments = (force);
indentation_increments = gradient(indentation.^beta);

indentation_velocity = medfilt1(indentation_increments ./ time_step,5);
cumsum_indentation_velocity = (indentation_velocity); % cumsum

relaxation_young_modulus = force_increments * (2 * n) / (time(end) - time(1)) / alpha;
relaxation_young_modulus = flip(relaxation_young_modulus ./ cumsum_indentation_velocity);

% Расчет релаксационного модуля Юнга на основе вейвлет-преобразования для
% индентора конической формы
% min_scale = 1;
% max_scale = 80;
% w0 = 10; % floor(indentation(end) / max_scale)
% n_scales = 2^9;
% 
% scale = zeros(1,n_scales);
% transform = zeros(2*n,n_scales);
% 
% for n_scale = 1:n_scales 
%   scale(n_scale) = exp(log((max_scale) / (min_scale)) * (n_scale - 1) / n_scales);
%   half_length = uint64(length(indentation) / 2); % 90
%   indentation_medium = indentation(half_length);
%   scaled_quantity_indentation = ((indentation - indentation_medium) ...
%    .* (indentation - indentation_medium)) ...
%    / (scale(n_scale) * scale(n_scale) * w0^2); % 
% 
%   indentation_norm = (indentation-indentation_medium) ./ (scale(n_scale) * w0);
%   Eexp_term = exp(-scaled_quantity_indentation ./ 2);
%   E1 = sqrt(scaled_quantity_indentation);
%   E2 = (scaled_quantity_indentation - 1) / (scale(n_scale) * scale(n_scale) * w0 * w0); % 
% 
%   gaussian_function_order_0 = Eexp_term / sum(Eexp_term); % 
% %   gaussian_function_order_1 = gaussian_function_order_0 .* E1; % TODO: не верно
%   gaussian_function_order_2 = gaussian_function_order_0 .* E2;
% 
%   force_with_flip = [flip(force); force]; % circshift(force,half_length)
% 
%   gaussian_function_fft = fft([zeros(half_length,1); gaussian_function_order_2; zeros(half_length,1)]);
%   force_fft = fft(force_with_flip); %  (half_length+1:3*half_length)
% 
%   transform(:,n_scale) = ifft(gaussian_function_fft .* force_fft);
% end
% 
% relaxation_young_modulus = max(((transform(1:end/2,:))),[],1);

end

