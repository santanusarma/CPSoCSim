% PLOT_POWER_MODEL plots the dynamic power vs speed, and the static power
% vs temperature for a specfied block on a core assigned to a specified
% thread. Optionally, it also plots the power profiles for the entire core
% corresponding to the specfied thread, and the total chip power. Note that
% this function assumes that leakage dependence on temperature is not
% ignored. 
function plot_power_model(block, thread, Power, Thread, ...
    plot_core_total, plot_chip_total, n_core)

%% Initialize figure constants
figure_handle = 1;
figure_file_base = 'results/tmp/figure';

%% Check if block and thread are within bounds
P_d = Thread.Power.Dynamic(:, thread);
n_units_per_core = length(P_d);
if block < 1 || block > n_units_per_core
    error('plot_power_model: block argument must be between 1 and %d', ...
        n_units_per_core); 
end
if thread < 1 || thread > Thread.n_thread
    error('plot_power_model: thread argument must be between 1 and %d', ...
        Thread.n_thread); 
end

%% Power of a particular thread and block

% Dynamic power at max. speed
p_d_max = Thread.Power.Dynamic(block, thread);

% Active leakage power as a function of temperature
C_a = Thread.Power.ActLeak.C(:, block, thread);
M_a = Thread.Power.ActLeak.M(:, block, thread);

% Standby leakage power as a function of temperature
C_s = Power.Mean.StdLeak.C(:, block);
M_s = Power.Mean.StdLeak.M(:, block); 

% Temperature and speed ranges
T = Power.Temperature_range;
I_split = Power.I_PWL_LDT;
s = linspace(0, 1, 10);

% Construct power vs speed and temperature
P_d = p_d_max * s;
P_a = evaluate_PWL_model(T, I_split, C_a, M_a);
P_s = evaluate_PWL_model(T, I_split, C_s, M_s);
P_total = p_d_max + P_a; % Total active power at s = 1

%% Plot block power
figure(figure_handle);
clf;
set(gcf, 'DefaultAxesFontName', 'AvantGarde');
set(gcf, 'DefaultLineLineWidth', 2);
set(gcf, 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontSize', 12);


subplot(221), plot(s, P_d);
xlabel('Normalized speed');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Dynamic power of block %d, thread %d', ...
    block, thread);
title(title_str);

subplot(222), plot(T, P_a);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Active leakage of block %d, thread %d', ...
    block, thread);
title(title_str);

subplot(223), plot(T, 1000*P_s);
xlabel('Die temperature (C)');
ylabel('Power (mW)');
grid on;
axis tight;
title_str = sprintf('Standby leakage of block %d, thread %d', ...
    block, thread);
title(title_str);

subplot(224), plot(T, P_total);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf(['Total active power of block %d, \n' ...
    'thread %d at max. speed'], block, thread);
title(title_str);

figure_file = sprintf('%s%d.pdf', figure_file_base, figure_handle);
print(gcf, '-dpdf', figure_file);
figure_handle = figure_handle + 1;

%% Compute total power of core on which thread resides
if (nargin < 5) || ~plot_core_total, return; end;

% Dynamic power at max. speed
p_d_max = sum(Thread.Power.Dynamic(:, thread));

% Active and standby leakage powers
P_a = zeros(length(T), 1);
P_s = zeros(length(T), 1);
for block = 1:n_units_per_core
    % Active leakage power as a function of temperature
    C_a = Thread.Power.ActLeak.C(:, block, thread);
    M_a = Thread.Power.ActLeak.M(:, block, thread);
    P_a = P_a + evaluate_PWL_model(T, I_split, C_a, M_a);

    % Standby leakage power as a function of temperature
    C_s = Power.Mean.StdLeak.C(:, block);
    M_s = Power.Mean.StdLeak.M(:, block); 
    P_s = P_s + evaluate_PWL_model(T, I_split, C_s, M_s);
end

% Construct power vs speed and temperature
P_d = p_d_max * s;
P_total = p_d_max + P_a; % Total active power at s = 1

%% Plot thread power
figure(figure_handle);
clf;
set(gcf, 'DefaultAxesFontName', 'AvantGarde');
set(gcf, 'DefaultLineLineWidth', 2);
set(gcf, 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontSize', 12);

subplot(221), plot(s, P_d);
xlabel('Normalized speed');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Dynamic power of thread %d', thread);
title(title_str);

subplot(222), plot(T, P_a);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Active leakage of thread %d', thread);
title(title_str);

subplot(223), plot(T, P_s);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Standby leakage of thread %d', thread);
title(title_str);

subplot(224), plot(T, P_total);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Total active power of thread %d at max. speed', ...
    thread);
title(title_str);

figure_file = sprintf('%s%d.pdf', figure_file_base, figure_handle);
print(gcf, '-dpdf', figure_file);
figure_handle = figure_handle + 1;

%% Compute total chip power
if (nargin < 6) || ~plot_chip_total, return; end;

% Dynamic power at max. speed
p_d_max = sum(sum(Thread.Power.Dynamic(:, :)));

% Active leakage
P_a = zeros(length(T), 1);
P_s = zeros(length(T), 1);
for thread = 1:Thread.n_thread
    for block = 1:n_units_per_core
        % Active leakage power as a function of temperature
        C_a = Thread.Power.ActLeak.C(:, block, thread);
        M_a = Thread.Power.ActLeak.M(:, block, thread);
        P_a = P_a + evaluate_PWL_model(T, I_split, C_a, M_a);
    end
end

% Standby leakage
for block = 1:n_units_per_core
    % Standby leakage power as a function of temperature
    C_s = Power.Mean.StdLeak.C(:, block);
    M_s = Power.Mean.StdLeak.M(:, block); 
    P_s = P_s + evaluate_PWL_model(T, I_split, C_s, M_s);
end
% Compute standby power for each core that is in standby, assuming all
% threads are active.
P_s = P_s * (n_core - Thread.n_thread);

% Construct power vs speed and temperature
P_d = p_d_max * s;
P_total = p_d_max + P_a + P_s; % Total power at s = 1

%% Plot total chip power
figure(figure_handle);
clf;
set(gcf, 'DefaultAxesFontName', 'AvantGarde');
set(gcf, 'DefaultLineLineWidth', 2);
set(gcf, 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontSize', 12);

subplot(221), plot(s, P_d);
xlabel('Normalized speed');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Dynamic power with %d threads, %d cores', ...
    Thread.n_thread, n_core);
title(title_str);

subplot(222), plot(T, P_a);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Active leakage with %d threads, %d cores', ...
    Thread.n_thread, n_core);
title(title_str);

subplot(223), plot(T, P_s);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf('Standby leakage with %d threads, %d cores', ...
    Thread.n_thread, n_core);
title(title_str);

subplot(224), plot(T, P_total);
xlabel('Die temperature (C)');
ylabel('Power (W)');
grid on;
axis tight;
title_str = sprintf(['Total power at full speed\nwith %d threads,' ...
    ' %d cores'], Thread.n_thread, n_core);
title(title_str);

figure_file = sprintf('%s%d.pdf', figure_file_base, figure_handle);
print(gcf, '-dpdf', figure_file);
