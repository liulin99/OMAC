function data = simulate_year_modes_full(output_dir, n_days)
%SIMULATE_YEAR_MODES_FULL Generate one continuous year of hourly glucose data.
%
% This driver uses the ODE models from:
%   Liu et al., "PID-controller enhanced artificial beta-cells"
%
% This file is self-contained. The physiological and PID ODE functions are
% included as local functions at the end of this file.
%
% Run:
%   data = simulate_year_modes_full();
%   data = simulate_year_modes_full(pwd, 30);
%
% Outputs:
%   glucose_simulation_hourly.csv
%   glucose_simulation_hourly.mat
%   glucose_simulation_hourly.png
%
% Time is measured in minutes inside the ODE model. Glucose is in mM.

if nargin < 1 || isempty(output_dir)
    output_dir = pwd;
end

if nargin < 2 || isempty(n_days)
    n_days = 365;
end

validateattributes( ...
    n_days, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
    mfilename, 'n_days');

if ~isfolder(output_dir)
    mkdir(output_dir);
end

rng(20260723, "twister");

minutes_per_hour = 60;
minutes_per_day = 24 * minutes_per_hour;
final_time = n_days * minutes_per_day;

% State order in the paper:
% [G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP]
glucose_index = 1;
stomach_glucose_index = 3;

%% Cached initial states copied from the comments in paper_code.m
X2_healthy = [ ...
    5.53723356890282; 25; 5.53723389333646; ...
    0.0423894127678293; 0.710387799226937; 0; ...
    12.5510163452660; 157.756021217255; ...
    0.0308742396866240; 0.103499263096041; ...
    0.585284413977824; 0.466763624715294; ...
    73.0063888613527; 3.17207489695228; ...
    0.267904086138662];

X2_implant = [ ...
    25.6264216894093; 25.6264216897280; ...
    25.6264216901862; 0.0424158564766771; ...
    0.00189884114381980; 0.00189655792181667; ...
    0.0335484296082557; 510.257201386281; ...
    0.0308742395239717; 0.103498772206185; ...
    0.586267230835652; 0.466763031906142; ...
    73.0130618152451; 3.17374557909211; ...
    0.268585211098882];

% Use a 16th state for the PID integral. In non-PID modes it is frozen.
state = [X2_healthy; 0];

%% Define the four biological modes
% The disease levels F_I = 0.5 and F_I = 0.1 are simulation choices.
% The beta-cell treatment is the original non-PID synthetic-cell model.
modes(1) = make_mode( ...
    "healthy", 1.0, 0, 0, 0, false);
modes(2) = make_mode( ...
    "moderate_untreated", 0.5, 0, 0, 0, false);
modes(3) = make_mode( ...
    "severe_untreated", 0.1, 0, 0, 0, false);
modes(4) = make_mode( ...
    "severe_beta_cell_treatment", ...
    0.1, 0.9137, 0.7047, 2e8, false);

% PID parameters used in the paper.
pid.Kp = 14;
pid.Ki = 11;
pid.Kd = 6;
pid.ATP_target = 3.2;

%% Simple one-year mode schedule
% Months 1-3: healthy
% Months 4-6: moderate disease, F_I = 0.5
% Months 7-9: severe disease, F_I = 0.1
% Months 10-12: severe disease with non-PID beta-cell treatment
schedule_boundaries = n_days * [0; 3; 6; 9; 12] / 12;

schedule_start_day = schedule_boundaries(1:end-1);
schedule_end_day   = schedule_boundaries(2:end);
schedule_mode_id   = [1; 2; 3; 4];

assert(schedule_start_day(1) == 0);
assert(abs(schedule_end_day(end) - n_days) < 10 * eps(max(1, n_days)));
assert(all(abs(schedule_start_day(2:end) - ...
    schedule_end_day(1:end-1)) < 10 * eps(max(1, n_days))));

%% Three reproducible oral-glucose inputs per day
% The original Stage 3 code changes state 3 (G_I) to represent the oral
% glucose challenge. We use the same mechanism at 08:00, 13:00, and 19:00.
meal_hours = [8, 13, 19];
meal_ranges_mM = [45, 65; 55, 80; 60, 90];

n_meals = n_days * numel(meal_hours);
meal_times = zeros(n_meals, 1);
meal_levels = zeros(n_meals, 1);
row = 0;
for day = 0:(n_days - 1)
    for meal = 1:numel(meal_hours)
        row = row + 1;
        meal_times(row) = ...
            day * minutes_per_day + meal_hours(meal) * minutes_per_hour;
        low = meal_ranges_mM(meal, 1);
        high = meal_ranges_mM(meal, 2);
        meal_levels(row) = low + (high - low) * rand;
    end
end

%% Event-driven continuous simulation
hour_times = (0:minutes_per_hour:final_time)';
mode_switch_times = schedule_start_day(2:end) * minutes_per_day;

% Safety filter: never simulate beyond the requested horizon.
mode_switch_times = mode_switch_times( ...
    mode_switch_times > 0 & mode_switch_times < final_time);

event_times = unique([0; meal_times; mode_switch_times; final_time]);
assert(all(event_times >= 0 & event_times <= final_time), ...
    'An event lies outside the requested simulation horizon.');

n_hours = numel(hour_times);
hourly_states = nan(n_hours, 16);
recorded = false(n_hours, 1);

ode_options = odeset( ...
    "NonNegative", 1:15, ...
    "RelTol", 1e-5, ...
    "AbsTol", 1e-8, ...
    "MaxStep", 10);

previous_mode_id = mode_at_time( ...
    0, schedule_start_day, schedule_end_day, schedule_mode_id, ...
    minutes_per_day);

% Runtime and remaining-time estimation
run_timer = tic;
next_progress_time = minutes_per_day;  % report once per simulated day

for event_index = 1:numel(event_times)
    t_start = event_times(event_index);
    current_mode_id = mode_at_time( ...
        t_start, schedule_start_day, schedule_end_day, ...
        schedule_mode_id, minutes_per_day);

    % At disease/treatment switches, blood and metabolic states continue.
    % When a new implant is introduced, initialize only its synthetic-cell
    % subsystem from the cached implant state.
    if current_mode_id ~= previous_mode_id
        entering_implant = ...
            modes(current_mode_id).N_S > 0 && ...
            modes(previous_mode_id).N_S == 0;
        if entering_implant
            implant_state_indices = [2, 4, 6, 9:15];
            state(implant_state_indices) = ...
                X2_implant(implant_state_indices);
            state(16) = 0;
        elseif modes(current_mode_id).use_pid && ...
                ~modes(previous_mode_id).use_pid
            state(16) = 0;
        end
    end

    % Apply a meal at the event time. max() prevents a meal event from
    % artificially lowering G_I if stomach glucose is already higher.
    meal_index = find(meal_times == t_start, 1);
    if ~isempty(meal_index)
        state(stomach_glucose_index) = max( ...
            state(stomach_glucose_index), meal_levels(meal_index));
    end

    % Record an hourly state after any instantaneous event at this time.
    hour_index = find(hour_times == t_start, 1);
    if ~isempty(hour_index)
        hourly_states(hour_index, :) = state';
        recorded(hour_index) = true;
    end

    if event_index == numel(event_times)
        break
    end

    t_end = event_times(event_index + 1);

    % Ask ode23 for any hourly samples strictly inside this event interval.
    % The endpoint is carried to the next event, where switches/meals are
    % applied before that hour is recorded.
    interior_hours = hour_times(hour_times > t_start & hour_times < t_end);
    t_eval = unique([t_start; interior_hours; t_end]);

    mode = modes(current_mode_id);
    rhs = @(t, x) mode_rhs(t, x, mode, pid);
    [t_segment, x_segment] = ode23(rhs, t_eval, state, ode_options);

    if ~isempty(interior_hours)
        [~, segment_rows] = ismember(interior_hours, t_segment);
        [~, output_rows] = ismember(interior_hours, hour_times);
        hourly_states(output_rows, :) = x_segment(segment_rows, :);
        recorded(output_rows) = true;
    end

    state = x_segment(end, :)';
    previous_mode_id = current_mode_id;
    
    % Display progress and estimated remaining runtime
    if t_end >= next_progress_time || t_end >= final_time
        elapsed_seconds = toc(run_timer);
        completed_fraction = t_end / final_time;
    
        estimated_total_seconds = ...
            elapsed_seconds / max(completed_fraction, eps);
    
        remaining_seconds = max( ...
            estimated_total_seconds - elapsed_seconds, 0);
    
        fprintf( ...
            "Simulated day %6.1f of %d | %5.1f%% complete | Elapsed: %6.2f hours | Estimated remaining: %6.2f hours\n", ...
            t_end / minutes_per_day, ...
            n_days, ...
            100 * completed_fraction, ...
            elapsed_seconds / 3600, ...
            remaining_seconds / 3600);
    
        % Schedule the next progress message
        while next_progress_time <= t_end
            next_progress_time = ...
                next_progress_time + minutes_per_day;
        end
    end

end

assert(all(recorded), "Not all hourly samples were recorded.");

%% Add temporally dependent observation noise at the hourly sampling rate
% stationary SD = 0.5 mM (about 9 mg/dL), AR coefficient = 0.7
rho = 0.7;
noise_sd_mM = 0.5;
innovation_sd_mM = noise_sd_mM * sqrt(1 - rho^2);
measurement_error = zeros(n_hours, 1);
measurement_error(1) = noise_sd_mM * randn;
for k = 2:n_hours
    measurement_error(k) = ...
        rho * measurement_error(k - 1) + innovation_sd_mM * randn;
end

glucose_true_mM = hourly_states(:, glucose_index);
glucose_observed_mM = max(glucose_true_mM + measurement_error, 0);
mM_to_mg_dL = 18.0182;

mode_id = zeros(n_hours, 1);
for k = 1:n_hours
    mode_id(k) = mode_at_time( ...
        hour_times(k), schedule_start_day, schedule_end_day, ...
        schedule_mode_id, minutes_per_day);
end
mode_name = strings(n_hours, 1);
for k = 1:n_hours
    mode_name(k) = modes(mode_id(k)).name;
end

%% Build hourly meal-input vector
meal_input_mM = zeros(n_hours, 1);

[is_meal_hour, meal_location] = ismember(hour_times, meal_times);

meal_input_mM(is_meal_hour) = ...
    meal_levels(meal_location(is_meal_hour));

%% Ensure every table variable is a column vector
time_hour = hour_times(:) / minutes_per_hour;
time_day = hour_times(:) / minutes_per_day;

mode_id = mode_id(:);
mode_name = mode_name(:);
meal_input_mM = meal_input_mM(:);
glucose_true_mM = glucose_true_mM(:);
glucose_observed_mM = glucose_observed_mM(:);
measurement_error = measurement_error(:);

glucose_true_mg_dL = ...
    glucose_true_mM * mM_to_mg_dL;

glucose_observed_mg_dL = ...
    glucose_observed_mM * mM_to_mg_dL;

%% Verify that all variables have the same number of rows
expected_rows = n_hours;

assert(numel(time_hour) == expected_rows);
assert(numel(time_day) == expected_rows);
assert(numel(mode_id) == expected_rows);
assert(numel(mode_name) == expected_rows);
assert(numel(meal_input_mM) == expected_rows);
assert(numel(glucose_true_mM) == expected_rows);
assert(numel(glucose_observed_mM) == expected_rows);
assert(numel(glucose_true_mg_dL) == expected_rows);
assert(numel(glucose_observed_mg_dL) == expected_rows);
assert(numel(measurement_error) == expected_rows);

%% Save raw simulation results before constructing the table
raw_backup_path = fullfile( ...
    output_dir, ...
    'glucose_simulation_raw_backup.mat');

save( ...
    raw_backup_path, ...
    'hour_times', ...
    'hourly_states', ...
    'mode_id', ...
    'mode_name', ...
    'meal_input_mM', ...
    'glucose_true_mM', ...
    'glucose_observed_mM', ...
    'measurement_error', ...
    'modes', ...
    'pid', ...
    'schedule_start_day', ...
    'schedule_end_day', ...
    'schedule_mode_id', ...
    'meal_times', ...
    'meal_levels', ...
    '-v7.3');

%% Construct output table
data = table( ...
    time_hour, ...
    time_day, ...
    mode_id, ...
    mode_name, ...
    meal_input_mM, ...
    glucose_true_mM, ...
    glucose_observed_mM, ...
    glucose_true_mg_dL, ...
    glucose_observed_mg_dL, ...
    measurement_error, ...
    'VariableNames', { ...
        'time_hour', ...
        'time_day', ...
        'true_mode_id', ...
        'true_mode_name', ...
        'meal_GI_input_mM', ...
        'glucose_true_mM', ...
        'glucose_observed_mM', ...
        'glucose_true_mg_dL', ...
        'glucose_observed_mg_dL', ...
        'measurement_error_mM'});

csv_path = fullfile(output_dir, "glucose_simulation_hourly.csv");
mat_path = fullfile(output_dir, "glucose_simulation_hourly.mat");
png_path = fullfile(output_dir, "glucose_simulation_hourly.png");

writetable(data, csv_path);
save(mat_path, ...
    "data", "hourly_states", "modes", "pid", ...
    "schedule_start_day", "schedule_end_day", "schedule_mode_id", ...
    "meal_times", "meal_levels");

%% Compact diagnostic plot
figure("Color", "w", "Position", [100, 100, 1400, 650]);
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
plot(data.time_day, data.glucose_true_mM, ...
    "k-", "LineWidth", 1.0);
hold on;
plot(data.time_day, data.glucose_observed_mM, ...
    "Color", [0.55, 0.75, 0.95], "LineWidth", 0.5);
hold off;
xlabel("Day");
ylabel("Glucose (mM)");
legend("True", "Observed", "Location", "best");
title(sprintf("Continuous %g-day glucose trajectory", n_days));
grid on;

nexttile;
stairs(data.time_day, data.true_mode_id, ...
    "LineWidth", 1.5, "Color", [0.65, 0.15, 0.15]);
xlabel("Day");
ylabel("True mode");
yticks(1:numel(modes));
yticklabels(string({modes.name}));
ylim([0.5, numel(modes) + 0.5]);
grid on;

exportgraphics(gcf, png_path, "Resolution", 160);
close(gcf);

fprintf('Saved %d hourly samples, from hour 0 through hour %.0f.\n', ...
    height(data), data.time_hour(end));
fprintf("CSV: %s\n", csv_path);
fprintf("MAT: %s\n", mat_path);
fprintf("PNG: %s\n", png_path);
end

function mode = make_mode(name, F_I, gamma_I, gamma_G, N_S, use_pid)
mode.name = name;
mode.F_I = F_I;
mode.gamma_I = gamma_I;
mode.gamma_G = gamma_G;
mode.N_S = N_S;
mode.use_pid = use_pid;
end

function id = mode_at_time( ...
        t_min, start_day, end_day, schedule_id, minutes_per_day)
day = t_min / minutes_per_day;
id_index = find(day >= start_day & day < end_day, 1, "last");
if isempty(id_index)
    % The final sample is exactly at the right endpoint, day 365.
    id_index = numel(schedule_id);
end
id = schedule_id(id_index);
end

function dx = mode_rhs(t, x, mode, pid)
if mode.use_pid
    dx = PID_model( ...
        t, x, mode.F_I, mode.gamma_I, mode.gamma_G, ...
        pid.Kp, pid.Ki, pid.Kd, pid.ATP_target, mode.N_S);
else
    dx15 = vivo_model( ...
        t, x(1:15), mode.F_I, ...
        mode.gamma_I, mode.gamma_G, mode.N_S);
    dx = [dx15; 0];
end
end
function dVarsdt = vivo_model(t, Vars, F_I, gamma_I, gamma_G, N_S)
    %%constants
    g_p = 0.1061; %mM/min glucose production; Table S10
    g_r = 4.1*10^(-3); %mM/min glucose effectiveness, insulin-independent glucose removal; Table S10
    k_si = 1.2*10^(-3); %L/ug/min, Insulin sensitivity;  Table S10
    % V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    % D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    % V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    % V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_iar = 0.0566; %/min, Insulin clearance from action compartment; Table S10
    k_iir = 0.0351; %/min, Insulin clearance from virtual compartment; Table S10
    %F_I = 0.01; %factor for diabetic insulin production; Table S10 !!!
    % k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    V_max1 = 10^(-13); %mmol/min/cell, maximum specific rate of cellular glucose uptake; Table S7
    K_m1 = 0.31; %mM, Michaelis-Menten constant for the updtake of glucose; Table S7
    F_SCI = 0.0229; %ug/U, scaling factor of insulin; Table S6
    k_tl = 5*10^(-8); %U/mM/min/cell, common translation constant; Table S9
    k_di = 1.1*10^(-3); %/min, Degradation constant insulin in the synthetic compartment; Table S9
    n_NFAT = 9; %no unit, Number of NFAT repeats for calcium dependent Insulin expression; Table S9
    n_S = 3.6; %no unit, plasmid copy number for SEAP or Insulin transfection; Table S9 and S11
    n_C = 3.6; %no unit, plasmid copy number for Cav1.3 transfection; Table S9 and S11
    V_max4 = 2.67*10^(-6); %mM/min, Max. calcium-dependent gene expression; Table S9
    d_4 = 1.6; %no unit, Hill coefficient for calcium-dependent gene expression; Table S9, not found in the paper, already sent email for inqury; unclear !!!
    K_m4 = 1.1*10^(-3); %mM, affinity for calcium-dependent gene expression; Table S9
    k_dm = 2.8*10^(-3); %/min, Degradation constant mRNA; Table S9
    r = 6.5; %mu_m, Radius of HEK-293 cell; Table S6
    V_c = 3.1415926*4*(r*10^(-5))^3/3; %L, volume of HEK-293 cell; Table S6
    F = 9.6485*10^(4); %C/mol, Faraday constant; Table S8
    R = 8.3144; %J/mol/K, ideal gas constant; Table S8
    T = 310; %K, temperature; Table S8
    C = (7.2/60)*10^(-3); %nF, Capacitance of HEK-293 cell (scaled for time units in min.); Table S8; unclear !!!
    alpha = 6*10^(-8)/(V_c*F); %mM/min/pA, Coefficient converting current in pA to concentration per unit time in mM/min; Table S8
    f_Ca = 10^(-3); %no unit, Fraction of cytosolic calcium that is free; Table S8
    g_Ca = 6.43*10^(-3); %nS, Conductance of calcium ions through background calcium channels; Table S8
    g_CaV = 0.8; %nS, Conductance of calcium ions through voltage activated calcium channels; Table S8
    I_CaP_max = 9.4; %pA, Maximum calcium pump current; Table S8
    g_K = 1.97*10^(-2); %nS, conductance of potassium ions through background potassium channels; Table S8
    g_KV = 1.24; %nS, conductance of poatssium ions through voltage-activated potassium channels; Table S8
    g_KATP = 1.96; %nS, conductance of poatssium ions through ATP-dependent potassium channels; Table S8
    d_K = 7; %no unit, hill coefficient of ATP-sensitive potassium channles; Table S8
    A_tot = 5; %mM, total concentration of ATP+ADP; Table S7
    K1 = 3; %mM, ADP binding constant, ATP-sensitive potassium channels; Table S8
    K2 = 2.21; %mM, ATP binding constant, ATP-sensitive potassium channels; Table S8
    g_NA = 4.37*10^(-2); %nS, Conductance of sodium ions through background sodium channels; Table S8
    I_Nak_max = 42.72; %pA, Maximum sodiumpotassium pump curren; Table S8
    K_K = 1; %/min, Rate constant of potassium uptake from the medium by NaK pump; Table S8
    K_Na = 11; %/min, Rate constant of sodium uptake from the cell by NaK pump; Table S8
    q = 1; %no unit; Autocatalytic stoichiometry; Table S7
    u = 0.12; %mM/min; Max. global rate of ATP consumption; Table S7
    k_delta = 2; %mM; Michaelis-Menten constant for ATP consumption; Table S7
    K_m2A = 1.61; %mM; Affinity constant for allosteric inhibition of ATP consumption by ATP; Table S7
    a = 1;%ATP hill coefficient; Table S7
    b = 1;%ATP hill coefficient; Table S7
    c = 1;%ATP hill coefficient; Table S7
    V_max2 = 10; %mM/min; Max. rate of intracellular glucose and ATP consumption; Table S7
    K_m2G = 1.95; %mM; Michaelis-Menten constant for glucose consumption; Table S7
    V_max3 = 9.96; %mM/min; Max. intermediate consumption rate; Table S7
    K_m3 = 1.83; %mM; Affinity for intermediary consumption inhibition by ATP; Table S7
    Ca_S = 1.8; %mM; Table S4 and page 9, roughly constant
    K_S = 5.3; %mM; Table S4 and page 9, roughly constant
    Na_S = 155; %mM; Table S4 and page 9, roughly constant
    % k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    % k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    % k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    %gamma_G = 0.7047;  %/min, Vascular exchange constant glucose; Table S10
    D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    %gamma_I = 0.9137; %/min; Vascular exchange constant insulin; Table S10
    %%end of constants
    dVarsdt = zeros(15,1);
    %define all 15 variables
    %Vars := [G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP] ;
    G = Vars(1);
    G_S = Vars(2);
    G_I = Vars(3);
    G_Si = Vars(4);
    I = Vars(5);
    I_S = Vars(6);
    I_a = Vars(7);
    I_i = Vars(8);
    M_S = Vars(9);
    Ca_si = Vars(10);
    K_Si = Vars(11);
    Na_Si = Vars(12);
    V = Vars(13);
    A = Vars(14);
    sixCP = Vars(15);
    %copy all the necessary equations
    %[G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP]
    %equation(1) is not needed in vivo; since mu will be zero in vivo(Table
    %S6),making mu_r(Ns) always zero in vivo
    % gamma_Gs = gamma_G*(G-G_S); %equation(37)
    gamma_GS = gamma_G*(G-G_S); %equation(37)
    R1_GS = G_S*V_max1/(G_S+K_m1); %equation(14)
    R_2_Gsi_A = (K_m2A^(2*b-a) * A^a * V_max2 * G_Si)/((K_m2A^(2*b)+A^(2*b))*(G_Si + K_m2G)); %Equation(15)
    dVarsdt(2) = -N_S*R1_GS + gamma_GS;%equation(2); dG_Sdt
    dVarsdt(4) = R1_GS/V_c - R_2_Gsi_A; %equation(3); dG_Sidt
    R_3_6CP_A = V_max3*(1-A/A_tot)*(K_m3^(2*c)*sixCP)/(A^(2*c+1)+K_m3^(2*c+1)); %Equation(16) !!!
    dVarsdt(15) = R_2_Gsi_A - R_3_6CP_A; %equation(4); dsixCPdt
    delta_A = u*A^2/(k_delta^2 + A^2); %Equation(18)
    dVarsdt(14) = -1*q*R_2_Gsi_A + (q+1)*R_3_6CP_A - delta_A; %equation(5); dAdt
    V_K_Ksi = 10^(3)*R*T*log(K_S/K_Si)/F; %equation(19)
    I_K_V_K_si = g_K*(V-V_K_Ksi); %equation(22)
    I_KV_V_K_si = g_KV*(V-V_K_Ksi) * ( 1/(1+exp((-15-V)/5.6)) ) * (1/(1+exp((43+V)/4.1))); %equation(23)
    I_KATP_V_KSi_A = g_KATP*(V-V_K_Ksi)*(1+ ((A_tot-A)/K1)^d_K )/(1+ ((A_tot-A)/K1)^d_K + (A/K2)^d_K); %equation(24)
    I_NaK_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28) unclear here without the V parameter, not defined in paper !!!
    I_NaK_V_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28)
    dVarsdt(11) = -1*alpha*(I_K_V_K_si + I_KV_V_K_si + I_KATP_V_KSi_A) + 2*alpha*I_NaK_Na_Si; %Equation(6) dK_Sidt
    V_Na_Nasi = 10^(3)*R*T*log(Na_S/Na_Si)/F; %equation(20)
    I_Na_V_NaSi = g_NA*(V-V_Na_Nasi); %equation(25)
    dVarsdt(12) = -1*alpha*I_Na_V_NaSi - 3*alpha*I_NaK_Na_Si; %Equation(7); dNa_sidt
    V_ca_Casi = 10^(3)*R*T*log(Ca_S/Ca_si)/F; %equation(21)
    I_Ca_V_Casi = g_Ca*(V-V_ca_Casi); %equation(26)
    I_Cav_V_Casi = n_C*g_CaV*(V - V_ca_Casi); %equation(27)
    I_CaP_V_Casi = I_CaP_max*Ca_si/(Ca_si+2*10^(-4));%equation(29)
    dVarsdt(10) = -0.5*alpha*f_Ca*(I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(8); dca_sidt
    dVarsdt(13) = (-1/C)*(I_K_V_K_si+I_KV_V_K_si+I_KATP_V_KSi_A+I_Na_V_NaSi+I_NaK_V_Na_Si+I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(9); dVdt
    R_4_CaSi = V_max4*(Ca_si)^d_4/(K_m4^d_4+(Ca_si)^d_4); %equation(17)
    dVarsdt(9) = n_NFAT*n_S*R_4_CaSi-(k_dm*M_S); %equation(10); dMsdt, mur(Ns) term is zero since mu is zero in vivo(Table S6)
    %equation(11); this is about the seap, not needed here
    gamma_IS = gamma_I*(I-I_S); %equation(36)
    dVarsdt(6) = F_SCI*k_tl*N_S*M_S - k_di*I_S + gamma_IS; %equation(12), dI_sdt
    %equation(13) is not needed since N_s will never change in vivo
    
    gamma_GI = D_G*(G-G_I); %equation(38)
    dVarsdt(1) = g_p - (g_r + k_si*I_a)*G - V_S*gamma_GS/V_B - V_B*gamma_GI/V_G; %equation(30), dGdt
    dVarsdt(3) = gamma_GI; %equation(31), dG_idt
    dVarsdt(7) = I - k_iar*I_a; %equation(32), dI_adt
    dVarsdt(8) = F_I*(G - k_iir*I_i); %equation(33), dI_idt, F_I value is unclear !!!
    n_HI = 8; %no unit
    k_MI = 10^(-2); %ng/ml in the original paper, which equals to µg/L in this paper
    cosi_I = I^n_HI/(I^n_HI + k_MI^n_HI); %unlabled, showed up on the very begining of page 14; unclear
    %Here is the guess:
    %https://www.cell.com/cms/10.1016/j.molcel.2014.06.007/attachment/992caabb-c4e2-4b4d-a3fb-9a0ddd8a1f60/mmc1.pdf
    %page 27
    %dVarsdt(5) = (F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I)*(1-I/20); %equation(34), dIdt; F_I value is unclear !!! add upper bount of 20mM for blood insulin
    dVarsdt(5) = (F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I)*(1-I/20);
    %equation(35), this is about the seap, not needed here
end

function dVarsdt = PID_model(t, Vars, F_I, gamma_I, gamma_G, p_gaining, i_gaining, d_gaining, steady_A, N_S)
    %%constants
    g_p = 0.1061; %mM/min glucose production; Table S10
    g_r = 4.1*10^(-3); %mM/min glucose effectiveness, insulin-independent glucose removal; Table S10
    k_si = 1.2*10^(-3); %L/ug/min, Insulin sensitivity;  Table S10
    V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_iar = 0.0566; %/min, Insulin clearance from action compartment; Table S10
    k_iir = 0.0351; %/min, Insulin clearance from virtual compartment; Table S10
    %F_I = 0.01; %factor for diabetic insulin production; Table S10 !!!
    k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    V_max1 = 10^(-13); %mmol/min/cell, maximum specific rate of cellular glucose uptake; Table S7
    K_m1 = 0.31; %mM, Michaelis-Menten constant for the updtake of glucose; Table S7
    F_SCI = 0.0229; %ug/U, scaling factor of insulin; Table S6
    k_tl = 5*10^(-8); %U/mM/min/cell, common translation constant; Table S9
    k_di = 1.1*10^(-3); %/min, Degradation constant insulin in the synthetic compartment; Table S9
    n_NFAT = 9; %no unit, Number of NFAT repeats for calcium dependent Insulin expression; Table S9
    n_S = 3.6; %no unit, plasmid copy number for SEAP or Insulin transfection; Table S9 and S11
    n_C = 3.6; %no unit, plasmid copy number for Cav1.3 transfection; Table S9 and S11
    V_max4 = 2.67*10^(-6); %mM/min, Max. calcium-dependent gene expression; Table S9
    d_4 = 1.6; %no unit, Hill coefficient for calcium-dependent gene expression; Table S9, not found in the paper, already sent email for inqury; unclear !!!
    K_m4 = 1.1*10^(-3); %mM, affinity for calcium-dependent gene expression; Table S9
    k_dm = 2.8*10^(-3); %/min, Degradation constant mRNA; Table S9
    r = 6.5; %mu_m, Radius of HEK-293 cell; Table S6
    V_c = 3.1415926*4*(r*10^(-5))^3/3; %L, volume of HEK-293 cell; Table S6
    F = 9.6485*10^(4); %C/mol, Faraday constant; Table S8
    R = 8.3144; %J/mol/K, ideal gas constant; Table S8
    T = 310; %K, temperature; Table S8
    C = (7.2/60)*10^(-3); %nF, Capacitance of HEK-293 cell (scaled for time units in min.); Table S8; unclear !!!
    alpha = 6*10^(-8)/(V_c*F); %mM/min/pA, Coefficient converting current in pA to concentration per unit time in mM/min; Table S8
    f_Ca = 10^(-3); %no unit, Fraction of cytosolic calcium that is free; Table S8
    g_Ca = 6.43*10^(-3); %nS, Conductance of calcium ions through background calcium channels; Table S8
    g_CaV = 0.8; %nS, Conductance of calcium ions through voltage activated calcium channels; Table S8
    I_CaP_max = 9.4; %pA, Maximum calcium pump current; Table S8
    g_K = 1.97*10^(-2); %nS, conductance of potassium ions through background potassium channels; Table S8
    g_KV = 1.24; %nS, conductance of poatssium ions through voltage-activated potassium channels; Table S8
    g_KATP = 1.96; %nS, conductance of poatssium ions through ATP-dependent potassium channels; Table S8
    d_K = 7; %no unit, hill coefficient of ATP-sensitive potassium channles; Table S8
    A_tot = 5; %mM, total concentration of ATP+ADP; Table S7
    K1 = 3; %mM, ADP binding constant, ATP-sensitive potassium channels; Table S8
    K2 = 2.21; %mM, ATP binding constant, ATP-sensitive potassium channels; Table S8
    g_NA = 4.37*10^(-2); %nS, Conductance of sodium ions through background sodium channels; Table S8
    I_Nak_max = 42.72; %pA, Maximum sodiumpotassium pump curren; Table S8
    K_K = 1; %/min, Rate constant of potassium uptake from the medium by NaK pump; Table S8
    K_Na = 11; %/min, Rate constant of sodium uptake from the cell by NaK pump; Table S8
    q = 1; %no unit; Autocatalytic stoichiometry; Table S7
    u = 0.12; %mM/min; Max. global rate of ATP consumption; Table S7
    k_delta = 2; %mM; Michaelis-Menten constant for ATP consumption; Table S7
    K_m2A = 1.61; %mM; Affinity constant for allosteric inhibition of ATP consumption by ATP; Table S7
    a = 1;%ATP hill coefficient; Table S7
    b = 1;%ATP hill coefficient; Table S7
    c = 1;%ATP hill coefficient; Table S7
    V_max2 = 10; %mM/min; Max. rate of intracellular glucose and ATP consumption; Table S7
    K_m2G = 1.95; %mM; Michaelis-Menten constant for glucose consumption; Table S7
    V_max3 = 9.96; %mM/min; Max. intermediate consumption rate; Table S7
    K_m3 = 1.83; %mM; Affinity for intermediary consumption inhibition by ATP; Table S7
    Ca_S = 1.8; %mM; Table S4 and page 9, roughly constant
    K_S = 5.3; %mM; Table S4 and page 9, roughly constant
    Na_S = 155; %mM; Table S4 and page 9, roughly constant
    k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    %gamma_G = 0.7047;  %/min, Vascular exchange constant glucose; Table S10
    D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    %gamma_I = 0.9137; %/min; Vascular exchange constant insulin; Table S10
    %%end of constants
    dVarsdt = zeros(16,1);
    %define all 15 variables
    %Vars := [G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP] ;
    G = Vars(1);
    G_S = Vars(2);
    G_I = Vars(3);
    G_Si = Vars(4);
    I = Vars(5);
    I_S = Vars(6);
    I_a = Vars(7);
    I_i = Vars(8);
    M_S = Vars(9);
    Ca_si = Vars(10);
    K_Si = Vars(11);
    Na_Si = Vars(12);
    V = Vars(13);
    A = Vars(14);
    sixCP = Vars(15);

    %copy all the necessary equations
    %[G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP]
    %equation(1) is not needed in vivo; since mu will be zero in vivo(Table
    %S6),making mu_r(Ns) always zero in vivo
    % gamma_Gs = gamma_G*(G-G_S); %equation(37)
    gamma_GS = gamma_G*(G-G_S); %equation(37)
    R1_GS = G_S*V_max1/(G_S+K_m1); %equation(14)
    R_2_Gsi_A = (K_m2A^(2*b-a) * A^a * V_max2 * G_Si)/((K_m2A^(2*b)+A^(2*b))*(G_Si + K_m2G)); %Equation(15)
    dVarsdt(2) = -N_S*R1_GS + gamma_GS;%equation(2); dG_Sdt
    dVarsdt(4) = R1_GS/V_c - R_2_Gsi_A; %equation(3); dG_Sidt
    %R_3_6CP_A = V_max3*(1-A/A_tot)*(K_m3^(2*c)*sixCP)/(A^(2*c)+K_m3^(2*c)); %Equation(16) !!!
    R_3_6CP_A = V_max3*(1-A/A_tot)*(K_m3^(2*c)*sixCP)/(A^(2*c+1)+K_m3^(2*c+1)); 
    dVarsdt(15) = R_2_Gsi_A - R_3_6CP_A; %equation(4); dsixCPdt
    delta_A = u*A^2/(k_delta^2 + A^2); %Equation(18)
    dVarsdt(14) = -1*q*R_2_Gsi_A + (q+1)*R_3_6CP_A - delta_A; %equation(5); dAdt
    V_K_Ksi = 10^(3)*R*T*log(K_S/K_Si)/F; %equation(19)
    I_K_V_K_si = g_K*(V-V_K_Ksi); %equation(22)
    I_KV_V_K_si = g_KV*(V-V_K_Ksi) * ( 1/(1+exp((-15-V)/5.6)) ) * (1/(1+exp((43+V)/4.1))); %equation(23)
    I_KATP_V_KSi_A = g_KATP*(V-V_K_Ksi)*(1+ ((A_tot-A)/K1)^d_K )/(1+ ((A_tot-A)/K1)^d_K + (A/K2)^d_K); %equation(24)
    I_NaK_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28) unclear here without the V parameter, not defined in paper !!!
    I_NaK_V_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28)
    dVarsdt(11) = -1*alpha*(I_K_V_K_si + I_KV_V_K_si + I_KATP_V_KSi_A) + 2*alpha*I_NaK_Na_Si; %Equation(6) dK_Sidt
    V_Na_Nasi = 10^(3)*R*T*log(Na_S/Na_Si)/F; %equation(20)
    I_Na_V_NaSi = g_NA*(V-V_Na_Nasi); %equation(25)
    dVarsdt(12) = -1*alpha*I_Na_V_NaSi - 3*alpha*I_NaK_Na_Si; %Equation(7); dNa_sidt
    V_ca_Casi = 10^(3)*R*T*log(Ca_S/Ca_si)/F; %equation(21)
    I_Ca_V_Casi = g_Ca*(V-V_ca_Casi); %equation(26)
    I_Cav_V_Casi = n_C*g_CaV*(V - V_ca_Casi); %equation(27)
    I_CaP_V_Casi = I_CaP_max*Ca_si/(Ca_si+2*10^(-4));%equation(29)
    dVarsdt(10) = -0.5*alpha*f_Ca*(I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(8); dca_sidt
    dVarsdt(13) = (-1/C)*(I_K_V_K_si+I_KV_V_K_si+I_KATP_V_KSi_A+I_Na_V_NaSi+I_NaK_V_Na_Si+I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(9); dVdt
    R_4_CaSi = V_max4*(Ca_si)^d_4/(K_m4^d_4+(Ca_si)^d_4); %equation(17)
    
    %dVarsdt(9) = n_NFAT*n_S*R_4_CaSi-(k_dm*M_S); %equation(10); dMsdt, mur(Ns) term is zero since mu is zero in vivo(Table S6)

    
    % steady_A = 3.08;%3.122852553601846; %get this value by 3days + 3 weeks + increase Gi and simulate for another 72 hours
    %tx1 = M_S-0.04;
    %cons = 10000;
    %dVarsdt(9) = 5*10^(-3)*(Vars(14)-steady_A) + 0*dVarsdt(14) + 0*(integral_error)-(k_dm*M_S)-(1+(exp(cons*tx1)-exp(-cons*tx1))/(exp(cons*tx1)+exp(-cons*tx1))); %mRNA is controlled by ATP  
    dVarsdt(16) = Vars(14)-steady_A; %the integration of ATP error
    %[Vars(9), Vars(14), (Vars(14)-steady_A), Vars(16), dVarsdt(14), (p_gaining*(Vars(14)-steady_A) + i_gaining*Vars(16) + d_gaining*(dVarsdt(14)) - (k_dm*M_S))*(1-Vars(9)/0.031)]
    dVarsdt(9) = (p_gaining*(Vars(14)-steady_A) + i_gaining*Vars(16) + d_gaining*(dVarsdt(14)) - (k_dm*M_S))*(1-Vars(9)/0.031);
    %dVarsdt(9) = p_gaining*(Vars(14)-steady_A) + i_gaining*Vars(16) + d_gaining*(dVarsdt(14)) - (k_dm*M_S);

    %equation(11); this is about the seap, not needed here
    gamma_IS = gamma_I*(I-I_S); %equation(36)
    dVarsdt(6) = F_SCI*k_tl*N_S*M_S - k_di*I_S + gamma_IS; %equation(12), dI_sdt
    %equation(13) is not needed since N_s will never change in vivo
    
    gamma_GI = D_G*(G-G_I); %equation(38)
    dVarsdt(1) = g_p - (g_r + k_si*I_a)*G - V_S*gamma_GS/V_B - V_B*gamma_GI/V_G; %equation(30), dGdt
    dVarsdt(3) = gamma_GI; %equation(31), dG_idt
    dVarsdt(7) = I - k_iar*I_a; %equation(32), dI_adt
    dVarsdt(8) = F_I*(G - k_iir*I_i); %equation(33), dI_idt, F_I value is unclear !!!
    n_HI = 8; %no unit
    k_MI = 10^(-2); %ng/ml in the original paper, which equals to µg/L in this paper
    cosi_I = I^n_HI/(I^n_HI + k_MI^n_HI); %unlabled, showed up on the very begining of page 14; unclear
    %Here is the guess:
    %https://www.cell.com/cms/10.1016/j.molcel.2014.06.007/attachment/992caabb-c4e2-4b4d-a3fb-9a0ddd8a1f60/mmc1.pdf
    %page 27
    tx2 = I-4.6383;
    %dVarsdt(5) = F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I - (1+(exp(cons*tx2)-exp(-cons*tx2))./(exp(cons*tx2)+exp(-cons*tx2)))*100 ; %equation(34), dIdt; F_I value is unclear !!!
    dVarsdt(5) = (F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I)*(1-I/20); %equation(34), dIdt; add upperbount of 20mM for blood insulin
    %equation(35), this is about the seap, not needed here
end
