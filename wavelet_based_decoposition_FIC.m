%%
n_scales = 2^10;

%% TEST EXAMPLE: a piecewise linear function
alpha = .1;

z = (1:2^12)';
z_contact = 800;

force = zeros(size(z));
force(1:z_contact) = alpha * z(1:z_contact);
force = flip(circshift(force,length(z) - z_contact));

scale_size = 10;
min_scale = 1;
max_scale = 80;
scales = exp(log( (max_scale) / (min_scale)) * ((1:n_scales)' - 1) / n_scales );
z_scales = scales * scale_size;

[wt_coefficients_0_order,wt_coefficients_1_order,wt_coefficients_2_order] = get_wavelet_transform_coefficients(force,z,scale_size,n_scales,max_scale,min_scale);
wt_coefficients_2_order(wt_coefficients_2_order < 0) = 0;

[maxima_line,indices] = max(wt_coefficients_2_order,[],1);
% n_edge_effects_points = 790;
% mask = ones(size(wt_coefficients_2_order));
% mask(1:n_edge_effects_points,:) = 0;
% maxima_lines = islocalmax(wt_coefficients_2_order,1) .* mask;
% indices = zeros(n_scales,1);
% maxima_line = zeros(n_scales,1);
% for n_scale = 1:n_scales
%   indices(n_scale) = find(maxima_lines(:,n_scale) == 1,1);
%   maxima_line(n_scale) = wt_coefficients_2_order(indices(n_scale),n_scale);
% end

line_width = 3;

% figure(1);
% plot(log10(z_scales),log10(maxima_line)); % (log10(z_scales) < 2)

figure(2);
plot(z_scales,max(wt_coefficients_2_order(:,:),[],1),'k','LineWidth',line_width) % n_edge_effects_points:end-n_edge_effects_points
ylabel('{\itR}, Па');xlabel('{\itZ}, нм');
set_figure;axis tight;

figure(3);
% mesh(wt_coefficients_2_order(:,1:end)'); % n_edge_effects_points:end-n_edge_effects_points
imagesc(wt_coefficients_2_order(1:end,1:end));
% mesh(log10(abs(wt_coefficients_2_order(:,1:end)')));
colormap('jet')
axis xy;
hold on; plot(indices,1:n_scales,'w')
xlim([500 1500])
xlabel('{\itb}, нм');ylabel('{\its}, нм');
c = colorbar;
c.Label.String = '{\itT}_{{\itg}^{(2)}}[{\itF}], Па';
set_figure;

figure(4);
plot(wt_coefficients_2_order(:,1),'k','LineWidth',line_width); % n_edge_effects_points:end-n_edge_effects_points
xlabel('{\itZ}, нм');ylabel('{\itT}_{{\itg}^{(2)}}[{\itF}], Па');
set_figure;axis tight;

% figure(5);
% mesh(wt_coefficients_1_order(:,1:end)');
% % imagesc(wt_coefficients_1_order(1:end,1:end)');
% colormap('jet')
% axis xy;

figure(6);
plot(wt_coefficients_1_order(:,1),'k','LineWidth',line_width);
xlabel('{\itZ}, нм');ylabel('{\itT}_{{\itg}^{(1)}}[{\itF}], нН/нм');
set_figure;axis tight;

% figure(7);
% mesh(wt_coefficients_0_order(1:end,1:end)');
% % imagesc(wt_coefficients_0_order(1:end,1:end)');
% colormap('jet')
% axis xy;

figure(8);hold on;
plot(force(:),'k','LineWidth',line_width);
plot(wt_coefficients_0_order(:,1),'k:','LineWidth',3*line_width);
xlabel('{\itZ}, нм');ylabel('{\itF}, нН');
legend({'{\itF}','{\itT}_{{\itg}^{(0)}}[{\itF}]'});
set_figure;axis tight;

%% ANALYZE OF THE KV-FRACTIONAL MODEL
E = 1000;
tau = 1e-1;
alpha = .75;
beta = .2;
lambda = 2;

time_step = .0014;
time_length = 10;
n_time_points = 1200;
time = linspace(time_step,time_length/2,n_time_points)';

F_loading = E * gamma(lambda + 1) * (time / time(end)).^lambda ...
  .* ( 1 / gamma(lambda + 1 - alpha) * (time / tau).^-alpha ...
      + 1 / gamma(lambda + 1 - beta) * (time / tau).^-beta );
relaxation_function_frac_model = E / gamma(1 - alpha) * (time / tau).^-alpha ...
  + E / gamma(1 - beta) * (time / tau).^-beta;
relaxation_function_frac_model = (relaxation_function_frac_model - relaxation_function_frac_model(end)) / (relaxation_function_frac_model(1) - relaxation_function_frac_model(end));

time = [time; linspace(time(end) + time_step,time_length,n_time_points)'];
F_loading = [flip(F_loading); zeros(size(F_loading))];

scale_size = 2 * time_step;
min_scale = 1;
max_scale = 4e2;
scales = exp(log( (max_scale) / (min_scale)) * ((1:n_scales)' - 1) / n_scales );
time_scales = scales * scale_size;

[wt_coefficients_0_order,wt_coefficients_1_order,wt_coefficients_2_order] = get_wavelet_transform_coefficients((F_loading),time,scale_size,n_scales,max_scale,min_scale);
wt_coefficients_2_order(wt_coefficients_2_order < 0) = 0;

n_edge_effects_points = 500;
mask = ones(size(wt_coefficients_2_order));
mask(1:n_edge_effects_points,:) = 0;

[maxima_line,indices] = max(wt_coefficients_2_order,[],1);
% maxima_lines = islocalmax(wt_coefficients_2_order,1) .* mask;
% indices = zeros(n_scales,1);
% maxima_line = zeros(n_scales,1);
% for n_scale = 1:n_scales
%   indices(n_scale) = find(maxima_lines(:,n_scale) == 1,1);
%   maxima_line(n_scale) = wt_coefficients_2_order(indices(n_scale),n_scale);
% end
maxima_line = (maxima_line - maxima_line(end)) / (maxima_line(1) - maxima_line(end));

line_width = 3;

% figure(1);
% plot(log10(z_scales(log10(z_scales) < 2)),log10(maxima_line(log10(z_scales) < 2)));

figure(2);hold on;
plot(time(1:length(relaxation_function_frac_model)),relaxation_function_frac_model,'k-','LineWidth',line_width);
% yyaxis right;
plot(time_scales,maxima_line,'k--','LineWidth',line_width);
xlim([time_scales(1) time_scales(end)]);
xlabel('{\itt}, с'); ylabel('{\itR}, МПа');
legend({'Аналитическое решение','Решение на основе вейвлет-преобразования'});
set_figure;

figure(3);
% mesh(wt_coefficients_2_order(1:end,1:end)'); % n_edge_effects_points:end-n_edge_effects_points
imagesc(wt_coefficients_2_order');
% mesh(log10(abs(wt_coefficients_2_order(:,1:end)')));
% imagesc(log10(abs(wt_coefficients_2_order(:,1:end)')));
colormap('jet')
axis xy;
hold on; plot(indices,1:n_scales,'w.-');
% xticklabels(round([time(end/5) time(2*end/5) time(3*end/5) time(4*end/5)]))
xlabel('{\itt}'); ylabel('{\its}');
c = colorbar;
c.Label.String = '{\itT}_{{\itg}^{(2)}}[{\itF}], MPa';
set_figure;

% figure(4);
% plot(time,wt_coefficients_2_order(:,1),'k','LineWidth',4); % n_edge_effects_points:end-n_edge_effects_points
% xlim([4.95 5.05]);
% xlabel('t, s'); ylabel(['T_{g^{(2)}}[F](t,s(1) = ' num2str(time_scales(1)) ' s), MPa']);
% set_figure;

% figure(5);
% mesh(wt_coefficients_1_order(:,1:end)');
% % imagesc(wt_coefficients_1_order(1:end,1:end)');
% colormap('jet')
% axis xy;
% xlabel('t'); ylabel('s');
% zlabel('T_{g^{(1)}}[F], nN/nm');
% set_figure;

% figure(6);
% plot(time,wt_coefficients_1_order(:,1));
% xlabel('t, s'); ylabel(['T_{g^{(1)}}[F](t,s(1) = ' num2str(time_scales(1)) ' s), nN/nm']);
% set_figure;

% figure(7);
% mesh(wt_coefficients_0_order(1:end,1:end)');
% % imagesc(wt_coefficients_0_order(1:end,1:end)');
% colormap('jet')
% axis xy;
% xlabel('t'); ylabel('s');
% zlabel('T_{g^{(0)}}[F], MPa');
% set_figure;

figure(8);hold on;
plot(time,F_loading(:),'k-','LineWidth',line_width);
plot(time,wt_coefficients_0_order(:,1),'k:','LineWidth',2*line_width);
xlabel('{\itt}, с'); ylabel('{\itF}, нН');
legend({'Дробная модель Кельвина-Фойгта',['{\itT}_{{\itg}^{(0)}}[{\itF}]({\itt},{\its}(1) = ' num2str(time_scales(1)) ' с)']});
set_figure;
%% ANALYZE AFM DATA (but with a spherical indenter)
path = 'D:\YandexDisk\Научная работа\Рабочий каталог\АСМ\data\2023-04-28\';
load([path 'cell08speed_.mat']);
load('force_ind_time_cell08speed_.mat');

n_force_curve = 1;

[~,index_max_indentation_time] = max(force(:,n_force_curve));
force_loading = force(1:index_max_indentation_time,n_force_curve);
z_displacement_loading = Data{1,1}(1:index_max_indentation_time)';

scale_size = 10;
min_scale = 1;
max_scale = 150;
n_signal = 2^floor(log2(index_max_indentation_time));
[wt_coefficients_0_order,wt_coefficients_2_order,scales] = get_wavelet_transform_coefficients(force_loading(1:n_signal),z_displacement_loading(1:n_signal),scale_size,n_scales,max_scale,min_scale);

figure(1);hold on;
% plot(z_displacement_loading(1:n_signal),force_loading(1:n_signal))
% plot(z_displacement_loading(1:n_signal),wt_coefficients_0_order(1:n_signal,1))
plot(scales,max((wt_coefficients_2_order(1:end,1:end)),[],1));

figure;
% mesh(wt_coefficients_0_order(1:end,1:end)');
% mesh(log10(abs(wt_coefficients_0_order(1:end,1:end)')));
% imagesc(wt_coefficients_2_order(1:end,1:end)');
% mesh(wt_coefficients_2_order(1:end,1:end)');
% imagesc(log10(abs(wt_coefficients_2_order(1:end,1:end)')));
mesh(log10(abs(wt_coefficients_2_order(1:end,1:end)')));
colormap('jet')
axis xy;

%% ANALYZE OF THE KV-FRACTIONAL MODEL
E = 1;
tau = 1e-1;
alpha = .75;
beta = .2;
lambda = 2;

time_contact = 800;
time = (1:time_contact)';

F_loading = E * gamma(lambda + 1) * (time / time(end)).^lambda ...
  .* ( 1 / gamma(lambda + 1 - alpha) * (time / tau).^-alpha ...
      + 1 / gamma(lambda + 1 - beta) * (time / tau).^-beta );
relaxation_function_frac_model = E / gamma(1 - alpha) * (time / tau).^-alpha ...
  + E / gamma(1 - beta) * (time / tau).^-beta;
% relaxation_function_frac_model = (relaxation_function_frac_model - relaxation_function_frac_model(end)) / (relaxation_function_frac_model(1) - relaxation_function_frac_model(end));

time = (1:2^12)';
F_loading = [flip(F_loading); zeros(size(time((time_contact+1):end)))];

scale_size = 10;
min_scale = 1;
max_scale = 80;
scales = exp(log( (max_scale) / (min_scale)) * ((1:n_scales)' - 1) / n_scales );
time_scales = scales * scale_size;

[wt_coefficients_0_order,wt_coefficients_1_order,wt_coefficients_2_order] = get_wavelet_transform_coefficients((F_loading),time,scale_size,n_scales,max_scale,min_scale);
wt_coefficients_2_order(wt_coefficients_2_order < 0) = 0;

relaxation_function_wt = flip(wt_coefficients_2_order(50:750,1));

square_velocity = relaxation_function_wt(1) / relaxation_function_frac_model(50);

line_width = 3;

% figure(1);
% plot(log10(z_scales(log10(z_scales) < 2)),log10(maxima_line(log10(z_scales) < 2)));

figure(2);hold on;
plot((time(50:750)),(relaxation_function_frac_model(50:750)),'k-','LineWidth',line_width*2);
plot((time(50:750)),relaxation_function_wt / square_velocity,'--','LineWidth',line_width);
% plot(log10(time(50:750)),log10(relaxation_function_frac_model(50:750)),'k-','LineWidth',line_width);
% plot(log10(time(50:750)),log10(relaxation_function_wt / square_velocity),'k--','LineWidth',line_width);
xlabel('{\itt}, с'); ylabel('{\itR}, Па');
legend({'Аналитическое решение','Решение на основе вейвлет-преобразования'},'Location','best');
set_figure;

figure(3);
% mesh(wt_coefficients_2_order(1:end,1:end)'); % n_edge_effects_points:end-n_edge_effects_points
imagesc(wt_coefficients_2_order');
% mesh(log10(abs(wt_coefficients_2_order(:,1:end)')));
% imagesc(log10(abs(wt_coefficients_2_order(:,1:end)')));
colormap('jet')
axis xy;
% hold on; plot(indices,1:n_scales,'w.-');
% xticklabels(round([time(end/5) time(2*end/5) time(3*end/5) time(4*end/5)]))
xlabel('{\itt}'); ylabel('{\its}');
c = colorbar;
c.Label.String = '{\itT}_{{\itg}^{(2)}}[{\itF}], Pa';
set_figure;

figure(4);
plot(time,wt_coefficients_2_order(:,1),'k','LineWidth',4); % n_edge_effects_points:end-n_edge_effects_points
% xlim([4.95 5.05]);
xlabel('t, s'); ylabel(['T_{g^{(2)}}[F](t,s(1) = ' num2str(time_scales(1)) ' s), Pa']);
set_figure;

% figure(5);
% mesh(wt_coefficients_1_order(:,1:end)');
% % imagesc(wt_coefficients_1_order(1:end,1:end)');
% colormap('jet')
% axis xy;
% xlabel('t'); ylabel('s');
% zlabel('T_{g^{(1)}}[F], nN/nm');
% set_figure;

% figure(6);
% plot(time,wt_coefficients_1_order(:,1));
% xlabel('t, s'); ylabel(['T_{g^{(1)}}[F](t,s(1) = ' num2str(time_scales(1)) ' s), nN/nm']);
% set_figure;

% figure(7);
% mesh(wt_coefficients_0_order(1:end,1:end)');
% % imagesc(wt_coefficients_0_order(1:end,1:end)');
% colormap('jet')
% axis xy;
% xlabel('t'); ylabel('s');
% zlabel('T_{g^{(0)}}[F], MPa');
% set_figure;

figure(8);hold on;
plot(time,F_loading(:),'k-','LineWidth',line_width);
plot(time,wt_coefficients_0_order(:,1),'k:','LineWidth',2*line_width);
xlabel('{\itt}, с'); ylabel('{\itF}, Па');
legend({'Дробная модель Кельвина-Фойгта',['{\itT}_{{\itg}^{(0)}}[{\itF}]({\itt},{\its}(1) = ' num2str(time_scales(1)) ' с)']});
set_figure;


