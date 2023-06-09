function [sigma] = get_sigma_with_relaxation_young_modulus(relaxation_young_modulus,deformation,time)
% Функция для расчета напряжений на основе релаксационного модуля Юнга.
%     Пример использования:
%     time_step = .013;
%     time = (time_step:time_step:10)';
%     time_loading = time(time <= 1);
% 
%     indentation = linspace(time_step,1,length(time_loading))';
%     
%     Fs = readmatrix('data.txt');
%     F = Fs(:,1);
%     F_loading = F(time <= 1);
%     
%     [relaxation_young_modulus] = get_relaxation_young_modulus([0; F_loading],[],[0; time_loading],3/2);
%     
%     [sigma] = get_sigma_with_relaxation_young_modulus(relaxation_young_modulus,[0;indentation],[0;time_loading]);
%     plot(time_loading,sigma)



% Проверка входных данных
if nargin < 1 || isempty(relaxation_young_modulus)
  error('Please pass relaxation_young_modulus values to the function.');
end

n = length(relaxation_young_modulus);
if sum(size(relaxation_young_modulus)) ~= (n + 1)
  error('The relaxation_young_modulus must be a vector (1D-array).');
end

if nargin < 3 || isempty(time)
  error('Please pass time values to the function.');
end

if nargin < 2 || isempty(deformation)
  error('Please pass deformation values to the function.');
end

% Расчет напряжений на основе метода трапеций при const шаге по времени
time_step = time(2) - time(1);

deformation_increments = diff(deformation);
deformation_velocity = deformation_increments ./ time_step;

sigma = flip(relaxation_young_modulus) .* deformation_velocity;
sigma = 2/3 * (time(end) - time(1)) / (2 * n) * sigma;

end

