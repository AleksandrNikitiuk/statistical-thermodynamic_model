function [...
  wavelet_transform_coefficients_0_order,...
  wavelet_transform_coefficients_1_order,...
  wavelet_transform_coefficients_2_order,scales]...
  = get_wavelet_transform_coefficients(...
  force,z_displacement,scale_size,n_scales,max_scale,min_scale)
% Функция для расчета коэффициентов вейвлет-преобразования силовой кривой.
%     Пример использования:
%     path = 'D:\YandexDisk\Научная работа\Рабочий каталог\АСМ\data\2023-04-28\';
%     load([path 'cell08speed_.mat']);
%     load('force_ind_time_cell08speed_.mat');
%     
%     n_force_curve = 1;
%     
%     [~,index_max_indentation_time] = max(force(:,n_force_curve));
%     force_loading = force(1:index_max_indentation_time,n_force_curve);
%     z_displacement_loading = Data{1,1}(1:index_max_indentation_time)';
%     
%     n_scale = 1;
%     n_signal = 2^floor(log2(index_max_indentation_time));
%     wt_coefficients_0_order = get_wavelet_transform_coefficients(force_loading(1:n_signal),z_displacement_loading(1:n_signal),[],n_scale);
%     
%     figure;hold on;
%     plot(z_displacement_loading(1:n_signal),force_loading(1:n_signal))
%     plot(z_displacement_loading(1:n_signal),wt_coefficients_0_order)

% Проверка входных данных
if nargin < 1 || isempty(force)
  error('Please pass force values to the function.');
end

n = length(force);
if sum(size(force)) ~= (n + 1)
  error('The force must be a vector (1D-array).');
end

if nargin < 6 || isempty(min_scale)
  min_scale = 1;
end

if nargin < 5 || isempty(max_scale)
  max_scale = 80;
end

if nargin < 4 || isempty(n_scales)
  n_scales = 2^8;
end

if nargin < 3 || isempty(scale_size)
  scale_size = 10;
end

if nargin < 2 || isempty(z_displacement)
  z_displacement = (1:n)';
end

% Инициализация переменных
scales = zeros(1,n_scales);
wavelet_transform_coefficients_0_order = zeros(3*n,n_scales);
wavelet_transform_coefficients_1_order = zeros(3*n,n_scales);
wavelet_transform_coefficients_2_order = zeros(3*n,n_scales);

half_length = uint64(length(z_displacement) / 2);
z_displacement_medium = z_displacement(half_length);
z_displacement_square = ((z_displacement - z_displacement_medium) ...
  .* (z_displacement - z_displacement_medium));

force_with_flip_values = [flip(force); force; flip(force)]; %  circshift(force,half_length + 1)
% force_with_flip_values = [linspace(force(1),0,length(force))'; force; linspace(force(end),0,length(force))']; %  circshift(force,half_length + 1)
force_ft = fft(force_with_flip_values); %  (half_length+1:3*half_length)

% Расчет коэффициентов вейвлет-преобразования силовой кривой
for n_scale = 1:n_scales
  
  scales(n_scale) = exp(log( (max_scale) / (min_scale)) * (n_scale - 1) / n_scales );

  scaled_z_displacement = (z_displacement - z_displacement_medium) / (scales(n_scale) * scale_size);
  scaled_z_displacement_squre = z_displacement_square ...
   / (scales(n_scale)^2 * scale_size^2);

  exp_term = exp(-scaled_z_displacement_squre ./ 2);
  exp_term_first_derivative_coefficients = -scaled_z_displacement / (scales(n_scale) * scale_size);
  exp_term_second_derivative_coefficients = (scaled_z_displacement_squre - 1) / (scales(n_scale)^2 * scale_size^2);

  gaussian_function_0_order = exp_term / sum(exp_term); %
  gaussian_function_1_order = gaussian_function_0_order .* exp_term_first_derivative_coefficients;
  gaussian_function_2_order = gaussian_function_0_order .* exp_term_second_derivative_coefficients;

  gaussian_function_0_order_ft = fft([zeros(2*half_length,1); gaussian_function_0_order; zeros(2*half_length,1)]);
  gaussian_function_1_order_ft = fft([zeros(2*half_length,1); gaussian_function_1_order; zeros(2*half_length,1)]);
  gaussian_function_2_order_ft = fft([zeros(2*half_length,1); gaussian_function_2_order; zeros(2*half_length,1)]);
  
  wavelet_transform_coefficients_0_order(:,n_scale) = ifft(gaussian_function_0_order_ft .* force_ft);
  wavelet_transform_coefficients_1_order(:,n_scale) = ifft(gaussian_function_1_order_ft .* force_ft);
  wavelet_transform_coefficients_2_order(:,n_scale) = ifft(gaussian_function_2_order_ft .* force_ft);

%   [~,f,as] = execute_fourier_analysis(gaussian_function_2_order,1);
%   [~,i] = max(as);
%   frequencies(n_scale) = f(i) / (scales(n_scale) * scale_size);
end

wavelet_transform_coefficients_0_order = circshift(wavelet_transform_coefficients_0_order,n/2+1,1);
wavelet_transform_coefficients_0_order = wavelet_transform_coefficients_0_order(1:n,:);
wavelet_transform_coefficients_1_order = circshift(wavelet_transform_coefficients_1_order,n/2+1,1);
wavelet_transform_coefficients_1_order = wavelet_transform_coefficients_1_order(1:n,:);
wavelet_transform_coefficients_2_order = circshift(wavelet_transform_coefficients_2_order,n/2+1,1);
wavelet_transform_coefficients_2_order = wavelet_transform_coefficients_2_order(1:n,:);

end

