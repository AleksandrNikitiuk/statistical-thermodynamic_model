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
%     [relaxation_young_modulus] = get_relaxation_young_modulus([0; F_loading],[],[0; time_loading],3/2);
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

force_increments = diff(force);
indentation_increments = diff(indentation.^beta);

indentation_velocity = indentation_increments ./ time_step;

relaxation_young_modulus = force_increments * (2 * n) / (time(end) - time(1)) / alpha;
relaxation_young_modulus = flip(relaxation_young_modulus ./ indentation_velocity);

end

