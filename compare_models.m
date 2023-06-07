n_force_curve = 1;

% time
time_step = .0013;
time = (time_step:time_step:10)';
time_loading = time(time <= 1);
time_dwell = time(time > 1);

Fs_fractional_model = readmatrix('data.txt') * 1000;
F_fractional_model = Fs_fractional_model(time > 1,n_force_curve);

F_slm = ...
  feval(fit(time_dwell,F_fractional_model,...
  'a+b*exp(-x/c)','StartPoint', [143, 196.8, 1.923]),time_dwell);
F_2_prony_elements = ...
  feval(fit(time_dwell,F_fractional_model,...
  'a+b*exp(-x/c)+d*exp(-x/e)'),time_dwell); % ,'StartPoint', [141.2, 150.3, 0.5543, 152, 2.348]
F_3_prony_elements = ...
  feval(fit(time_dwell,F_fractional_model,...
  'a+b*exp(-x/c)+d*exp(-x/e)+f*exp(-x/g)'),time_dwell); % ,'StartPoint', [141.2, 150.3, 0.5543, 152, 2.348, 1157, 0.3522]

% visualiztion
figure(1);hold on;
plot(time_dwell,F_fractional_model,'k','LineWidth',2);
plot(time_dwell,F_slm,'b','LineWidth',2);
plot(time_dwell,F_2_prony_elements,'g','LineWidth',2);
% plot(time_dwell,F_3_prony_elements,'r','LineWidth',2);
xlabel('time, s');
ylabel('{\itF}, N');
set_figure;

figure(2);hold on;
plot(time_dwell,F_fractional_model,'k','LineWidth',2);
plot(time_dwell,F_slm,'b','LineWidth',2);
% plot(time_dwell,F_2_prony_elements,'g','LineWidth',2);
plot(time_dwell,F_3_prony_elements,'r','LineWidth',2);
xlabel('time, s');
ylabel('{\itF}, N');
set_figure;

