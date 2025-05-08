clc
close all
clear

%% Question 1

mass = 1.8;  % kg
diam_prop = 1.2;  % m
rpm = 2800;
omega = rpm * 2 *pi / 60;
peak_load = 510;  % W
req_load = 360;  % W

number_of_blades = 2;

rho = 20e-3;  % kg/m^3
grav_acc = 3.73;  % m/s^2
Cd_0 = 0.02;

tip_speed_ing = omega*(diam_prop/2);

T = mass * grav_acc / 2;
P_ideal = sqrt(T^3) / sqrt(2*rho*pi*(diam_prop/2)^2);  % for 1 rotor

data = load('chord.txt');

r_div_R = data(:,1);
c_div_R = data(:,2);

needed_r_div_R = r_div_R(17);
needed_c_div_R = c_div_R(17);

finally_c = needed_c_div_R * diam_prop/2;
finally_r = needed_r_div_R * diam_prop/2;

P_0 = 1/8 * rho * finally_c * number_of_blades * Cd_0 * omega^3 * (diam_prop/2)^4;

gamma = 1.15;

P_tot = 2*(gamma * P_ideal + P_0); % Total Power


%% Question 2

nr_prop = [2, 4];
R_new = linspace(0.3, 1.0, 71);
new_batt_mass = 0.5;  % kg
old_batt_capaicty = 10*3600;  % J
new_batt_capaicty = 20*3600;  % J

n_blades = [2, 3, 4];

m_no_fuselage_interm = zeros(length(R_new), 2, 3);
m_fuselage_interm = zeros(length(R_new), 2, 3);
m_total_interm = zeros(length(R_new), 2, 3);

m_no_fuselage = zeros(length(R_new), 2, 3);
m_fuselage = zeros(length(R_new), 2, 3);
m_total = zeros(length(R_new), 2, 3);
m_propeller_q2 = zeros(length(R_new), 2, 3);
m_control_q2 = zeros(length(R_new), 2, 3);
m_computer_q2 = zeros(length(R_new), 2, 3);

omega_2 = zeros(length(R_new), 2, 3);

P_new = zeros(length(R_new), 2, 3);
T_new = zeros(length(R_new), 2, 3);

for k = 1:length(n_blades)
    for j = 1:length(nr_prop)
        for i = 1:length(R_new)

            % Initialize diff to enter the while loop
            diff = 20;
            eps = 1e-6;
            omega_old = omega;
            m_old = mass;
            P_old = P_tot;
            while diff > eps
                % First mass estimation with given ingenuity power (P_tot)
                [m_propeller_q2(i,j,k), m_control_q2(i,j,k), m_computer_q2(i,j,k), m_no_fuselage_interm(i, j, k), m_fuselage_interm(i, j, k), m_total(i,j,k)] = ...
                    calc_mass(n_blades(k), R_new(i), P_old, 2*P_tot, nr_prop(j), new_batt_mass, 2);
                [T_per_rotor, P_per_rotor] = calculate_Power_Thrust(m_total(i,j,k), grav_acc, nr_prop(j), rho, R_new(i), needed_c_div_R, n_blades(k), Cd_0, omega_old);
                P_new(i,j,k) = P_per_rotor * nr_prop(j);
                T_new(i,j,k) = T_per_rotor * nr_prop(j);
                omega_2(i,j,k) = tip_speed_ing / R_new(i);
                diff = abs(m_old - m_total(i,j,k));
                if diff > eps
                    m_old = m_total(i,j,k);
                    omega_old = omega_2(i,j,k);
                    P_old = P_per_rotor;
                end
            end
        end
    end
end


% Calculations
% For bicopter (nr_prop(1)) - 2 Blades shows optimal performance
[P_bi_min, idx_bi] = min(P_new(:, 1, 1));  % j = 1, k = 1
R_opt_bi = R_new(idx_bi);
P_opt_bi = P_new(idx_bi, 1, 1);
omega_opt_bi = omega_2(idx_bi, 1, 1);
m_total_opt_bi = m_total(idx_bi, 1, 1);  % j = 1, k = 1
flight_time_bi = (new_batt_capaicty/(P_opt_bi))/60; % In minutes

% For quadcopter (nr_prop(2)) - 2 Blades shows optimal performance
[P_quad_min, idx_quad] = min(P_new(:, 2, 1));  % j = 2, k = 1
R_opt_quad = R_new(idx_quad);
P_opt_quad = P_quad_min;
T_opt_quad = T_new(idx_quad, 2, 1);
omega_opt_quad_q2 = omega_2(idx_quad, 2, 1);
m_total_opt_quad = m_total(idx_quad, 2, 1);
flight_time_quad = (new_batt_capaicty/(P_opt_quad))/60; % In minutes


% Mass components (in kg)
m_control_final_2 = m_control_q2(idx_quad, 2, 1);
m_fuselage_final_2 = m_fuselage_interm(idx_quad, 2, 1);
m_battery_final_2 = 0.5; % Given
m_propellor_final_2 = m_propeller_q2(idx_quad, 2, 1);
m_computer_final_2 = m_computer_q2(idx_quad, 2, 1);

% Mass data and labels
mass_data = [m_control_final_2, m_fuselage_final_2, m_battery_final_2, ...
             m_propellor_final_2, m_computer_final_2];

labels = {'Control + Motors', 'Fuselage', 'Battery', 'Propellers', 'Computer & Other'};

% Plot pie chart
figure;
pie(mass_data, labels);
title('Mass Distribution of Quadcopter Components');
save_plot('Mass distribution')


% Print required output
fprintf('\n========= BICOPTER CONFIGURATION =========\n');
fprintf('Optimal Radius (R_opt_bi): %.3f m\n', R_opt_bi);
fprintf('Minimum Power (P_opt_bi): %.3f W\n', P_opt_bi);
fprintf('Optimal Angular Velocity (omega_opt_bi): %.3f rad/s\n', omega_opt_bi);
fprintf('Total Mass (m_total_opt_bi): %.3f kg\n', m_total_opt_bi);
fprintf('Estimated Flight Time: %.2f minutes\n', flight_time_bi);

fprintf('\n========= QUADCOPTER CONFIGURATION =========\n');
fprintf('Optimal Radius (R_opt_quad): %.3f m\n', R_opt_quad);
fprintf('Minimum Power (P_opt_quad): %.3f W\n', P_opt_quad);
fprintf('Optimal Angular Velocity (omega_opt_quad): %.3f rad/s\n', omega_opt_quad_q2);
fprintf('Total Mass (m_total_opt_quad): %.3f kg\n', m_total_opt_quad);
fprintf('Estimated Flight Time: %.2f minutes\n', flight_time_quad);


% Plot for Bicopter
figure;
hold on;
for i = 1 : 3
    plot(R_new, P_new(:, 1, i), 'LineWidth', 1.5);
end

% Mark the optimal point for bicopter (2 blades, j=1, k=1)
plot(R_opt_bi, P_opt_bi, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Radius (m)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
grid on
yl = ylim;  % Get current y-axis limits
ylim([P_opt_bi - 100, yl(2)]);

title('Bicopter')
legend({'2 blades', '3 blades', '4 blades', 'Optimal (2 blades)'}, 'Location', 'best');
hold off;
save_plot('Bicopter_Q2');

% Plot for Quadcopter
figure;
hold on;
for i = 1 : 3
    plot(R_new, P_new(:, 2, i), 'LineWidth', 1.5);
end

% Mark the optimal point for quadcopter (2 blades, j=2, k=1)
plot(R_opt_quad, P_opt_quad, 'go', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Radius (m)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
grid on
% Find the max power across all blade configurations for j=2
P_max = max(P_new(:, 2, :), [], 'all');  % Safe and robust
ylim([P_opt_quad - 50, P_max + 50]);

title('Quadcopter')
legend({'2 blades', '3 blades', '4 blades', 'Optimal (2 blades)'}, 'Location', 'best');
hold off;
save_plot('Quadcopter_Q2');

%% Task 3

each_battery_mass = 0.047; % kg
each_battery_capacity = (10/6)*3600; % J
n_add_battery = linspace(0,300,301);
n_battery = 6 + n_add_battery;

m_no_fuselage_3 = zeros(length(n_battery),1);
m_fuselage_3 = zeros(length(n_battery),1);
m_total_3 = zeros(length(n_battery),1);
Power_3 = zeros(length(n_battery),1);
Thrust_3 = zeros(length(n_battery),1);
flight_time_3 = zeros(length(n_battery),1);


for i = 1:length(n_battery)
    diff_3 = 20;
    eps = 1e-6;
    m_old_3 = m_total_opt_quad;
    P_old_3 = P_opt_quad/4;
    while diff_3 > eps
        m_battery = n_battery(i)*each_battery_mass;
        cap_battery = n_battery(i)*each_battery_capacity;
        [~, ~, ~, m_no_fuselage_3, m_fuselage_3, m_total_3(i)] = calc_mass(2, R_opt_quad, P_old_3, 2*P_tot, 4, m_battery, 0);
        [Thrust_3_per_rotor, Power_3_per_rotor] = calculate_Power_Thrust(m_total_3(i), grav_acc, 4, rho, R_opt_quad, needed_c_div_R, 2, Cd_0, omega_opt_quad_q2);
        Thrust_3(i) = 4*Thrust_3_per_rotor; % Multiplying by no. of propellors
        Power_3(i) = 4* Power_3_per_rotor; % Multiplying by no. of propellors
        diff_3 = abs(m_old_3 - m_total_3(i));
        if diff_3 > eps
            m_old_3 = m_total_3(i);
            P_old_3 = Power_3(i)/4;
        end
        flight_time_3(i) = cap_battery/(Power_3(i));
    end
end

% Find the highest flight time and its corresponding n_battery
[max_flight_time_quad, idx] = max(flight_time_3);  % idx gives the index of the max flight time
max_add_battery_num_quad = n_add_battery(idx);  % Get the corresponding n_battery
max_pow_batt = max_add_battery_num_quad + 6;

% Find the flight time with the limited amount of battery to be carried
limit_batt = 42;
limit_add_batt = 42 - 6;
idx_limit = find(n_add_battery == limit_add_batt);
flight_time_limit_batt = flight_time_3(idx_limit);
Thrust_opt_quad_battery = Thrust_3(idx_limit);
Power_opt_quad_battery = Power_3(idx_limit);
mass_quad_with_battery = m_total_3(idx_limit);


% Display the result
fprintf('\n========= Quadcopter Flight Time for Optimal no. of %d Batteries =========\n', max_pow_batt);
fprintf('Flight Time: %.2f minutes\n', max_flight_time_quad / 60);
fprintf('Total Batteries Used: %.0f\n', max_pow_batt);

% Display the result
fprintf('\n========= Quadcopter Flight Time for Limit of %d Batteries =========\n', limit_batt);
fprintf('Flight Time: %.2f minutes\n', flight_time_limit_batt / 60);
fprintf('Total Batteries Used: %.0f\n', limit_batt);


figure;
plot(n_add_battery, flight_time_3/60, 'b-', 'LineWidth', 1.5);
hold on;

% Vertical line at max flight time
xline(max_add_battery_num_quad, 'k--', 'LineWidth', 1.5);

% Vertical line at battery limit
xline(limit_batt - 6, 'r--', 'LineWidth', 1.5);

xlabel('No. of Additional Batteries');
ylabel('Flight Time (min)');
title('Flight Time vs. No. of Additional Batteries - Quadcopter');
grid on;
legend({'Flight time', 'Max flight time', 'Limitted (2 kg) flight time'}, 'Location', 'southeast');
hold off;
save_plot('Power_battery_Q3');


%% Task 4

% Mars parameters
rho_mars = 0.019;           % kg/m^3
visc_mars = 9.82e-6;        % Pa·s

% Earth parameters
rho_earth = 1.225;          % kg/m^3
visc_earth = 1.81e-5;       % Pa·s

% Rotor parameters
R = R_opt_quad;             % Rotor radius
r = 0.75 * R;               % 75% span location
w = 293.215314335;          % Rotational speed in rad/s
V = r * w;                  % Local blade speed
c = 0.1209;                 % Chord length

% Reynolds numbers
Re_mars = rho_mars * V * c / visc_mars;
Re_earth = rho_earth * V * c / visc_earth;

% Ratio
Re_ratio = Re_earth / Re_mars;

% Display results
disp(['Reynolds number on Mars: ', num2str(Re_mars)])
disp(['Reynolds number on Earth: ', num2str(Re_earth)])
disp(['Reynolds ratio (Earth/Mars): ', num2str(Re_ratio)])

% Lift and drag coefficients of E63 airfoil

% Open the file
fileID = fopen('E63_aero.txt', 'r');

% Skip the first line
fgetl(fileID);  % This reads the first line and discards it

% Read the rest of the data
data_E63 = textscan(fileID, '%f %f %f', 'Delimiter', '\t');  % Adjust the format based on your data

aoa_E63_cell = data_E63(:,1); % in deg
cl_E63_cell = data_E63(:,2);
cd_E63_cell = data_E63(:,3);

aoa_E63 = aoa_E63_cell{1};
cl_E63 = cl_E63_cell{1};
cd_E63 = cd_E63_cell{1};

% Assuming aoa_E63, cl_E63, and cd_E63 are vectors containing the data
target_cl = 1.2204;
target_cd = 0.02841;

% Find the index where cl_E63 is closest to target_cl and cd_E63 is closest to target_cd
[~, idx_cl] = min(abs(cl_E63 - target_cl));
[~, idx_cd] = min(abs(cd_E63 - target_cd));

% You can check if the indices match or if they are the same point
if idx_cl == idx_cd
    aoa_d_E63 = aoa_E63(idx_cl);
else
    % If they do not match, you can calculate a weighted average or other method to find the aoa
    % Here, I'm simply averaging the two aoa values from both indices
    aoa_d_E63 = mean([aoa_E63(idx_cl), aoa_E63(idx_cd)]);
end

disp(['The angle of attack at Cl = ', num2str(target_cl), ' and Cd = ', num2str(target_cd), ' is: ', num2str(aoa_d_E63)]);




% Enhanced Plot: CL vs CD
figure;
plot(cd_E63, cl_E63, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0 0.4470 0.7410]);
grid on;
xlabel('Drag Coefficient (C_D)', 'FontSize', 12);
ylabel('Lift Coefficient (C_L)', 'FontSize', 12);
title('Lift vs Drag Coefficient (E63 Airfoil)', 'FontSize', 14);
legend('E63 Airfoil', 'Location', 'best');
set(gca, 'FontSize', 11);
save_plot('cl_cd_Q4');

% Enhanced Plot: CL vs AOA
figure;
plot(aoa_E63, cl_E63, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.8500 0.3250 0.0980]);
grid on;
xlabel('Angle of Attack (°)', 'FontSize', 12);
ylabel('Lift Coefficient (C_L)', 'FontSize', 12);
title('Lift Coefficient vs Angle of Attack (E63 Airfoil)', 'FontSize', 14);
legend('E63 Airfoil', 'Location', 'best');
set(gca, 'FontSize', 11);
save_plot('vl_aoa_Q4');


%% Task 5

% Varying chord and constant twist

c_tip = linspace(0.01,0.02,5);
Area_opt_quad = pi*(R_opt_quad^2);
y = linspace(0.01,0.98,100);
chord_E63 = zeros(length(y),length(c_tip));
radial_pos_E63 = zeros(length(y),length(c_tip));
theta_E63 = zeros(length(y),length(c_tip));
dP_E63 = zeros(length(y),length(c_tip));
dT_E63 = zeros(length(y),length(c_tip));
dCT_E63 = zeros(length(y),length(c_tip));
cl_q5 = zeros(length(y),length(c_tip));
cd_q5 = zeros(length(y),length(c_tip));
aoa = zeros(length(y),length(c_tip));
T_total_Q5 = zeros(length(c_tip),1);
P_total_Q5 = zeros(length(c_tip),1);
CT_quad = T_opt_quad/(rho_mars *(Area_opt_quad)*((omega_opt_quad_q2*R_opt_quad)^2));
for j = 1:length(c_tip)
    for i = 1:length(y)
        
        radial_pos_E63(i,j) = y(i)*R_opt_quad;
        chord_E63(i,j) = chord_airfoil(c_tip(j), radial_pos_E63(i,j), R_opt_quad);
        theta_E63(i,j) = twist_distribution(radial_pos_E63(i,j), R_opt_quad, 18, 0, 2);
        [dT_E63(i,j), dP_E63(i,j), dCT_E63(i,j), aoa(i,j), cl_q5(i,j), cd_q5(i,j)]  = calculate_BEM_final(omega_opt_quad_q2, radial_pos_E63(i,j), theta_E63(i,j), aoa_E63, cl_E63, cd_E63, 2, rho_mars, chord_E63(i,j), R_opt_quad);
    end
    T_total_Q5(j) = 4*trapz(radial_pos_E63(:,j), dT_E63(:,j));   % in Newtons - Multiplied with no. of propellors
    P_total_Q5(j) = 4*trapz(radial_pos_E63(:,j), dP_E63(:,j));   % in Watts - - Multiplied with no. of propellors
end
% c_E63_75 = 0.0510;
% Enhanced plot for chord distribution
figure;
hold on;
for j = 1:length(c_tip)
    plot(y, chord_E63(:,j), 'LineWidth', 2);
end
xlabel('Normalized Radius (y)', 'FontSize', 12);
ylabel('Chord length [m]', 'FontSize', 12);
title('Chord Distribution Along Blade Span', 'FontSize', 14);
legend(arrayfun(@(x) sprintf('c_{tip} = %.3fm', x), c_tip, 'UniformOutput', false), 'Location', 'best');
grid on;
xlim([min(y), max(y)]);
set(gca, 'FontSize', 12);
save_plot('Chord_distribution_Q5');

figure;
hold on;
for j = 1:length(c_tip)
    plot(y, rad2deg(aoa(:,j)), 'LineWidth', 2);
end
xlabel('Normalized Radius (-)', 'FontSize', 12);
ylabel('Angle of attack (deg)', 'FontSize', 12);
title('AoA Distribution Along Blade Span', 'FontSize', 14);
legend(arrayfun(@(x) sprintf('c_{tip} = %.3fm', x), c_tip, 'UniformOutput', false), 'Location', 'best');
grid on;
xlim([min(y), max(y)]);
set(gca, 'FontSize', 12);
save_plot('AoA_distribution_Q5');



% Enhanced plot for twist distribution
figure;
plot(y, theta_E63, 'r-', 'LineWidth', 2);
xlabel('Normalized Radius (-)', 'FontSize', 12);
ylabel('Twist angle \theta [deg]', 'FontSize', 12);
title('Twist Distribution Along Blade Span', 'FontSize', 14);
grid on;
xlim([min(y), max(y)]);
set(gca, 'FontSize', 12);
save_plot('Twist_distribution_Q5');



% Plot dP vs y
figure;
hold on;
for j = 1:length(c_tip)
    plot(y, dP_E63(:,j), 'LineWidth', 2);
end
xlabel('Normalized Radius (-)');
ylabel('Differential Power (dP)');
title('Differential Power Distribution along Blade Span');
legend(arrayfun(@(x) sprintf('c_{tip} = %.3fm', x), c_tip, 'UniformOutput', false), 'Location', 'best');
grid on;
xlim([min(y), max(y)]);
ylim([0, max(dP_E63, [], 'all') * 1.1]);
set(gca, 'FontSize', 12);
save_plot('dP_distribution_Q5');


% Plot dT vs y
figure;
hold on;
for j = 1:length(c_tip)
    plot(y, dT_E63(:,j), 'LineWidth', 2);
end
xlabel('Normalized Radius (y)');
ylabel('Differential Thrust (dT)');
title('Differential Thrust Distribution along Blade Span');
legend(arrayfun(@(x) sprintf('c_{tip} = %.3fm', x), c_tip, 'UniformOutput', false), 'Location', 'best');
grid on;
xlim([min(y), max(y)]);
ylim([0, max(dT_E63, [], 'all') * 1.1]);
set(gca, 'FontSize', 12);
save_plot('dT_distribution_Q5');


% Plot dCT vs y
figure;
hold on;
for j = 1:length(c_tip)
    plot(y, dCT_E63(:,j), 'LineWidth', 2);
end
xlabel('Normalized Radius (y)');
ylabel('Differential Thrust Coefficient (dC_T)');
title('Differential Thrust Coefficient Distribution along Blade Span');
legend(arrayfun(@(x) sprintf('c_{tip} = %.3fm', x), c_tip, 'UniformOutput', false), 'Location', 'best');
grid on;
xlim([min(y), max(y)]);
ylim([0, max(dCT_E63, [], 'all') * 1.1]);
set(gca, 'FontSize', 12);
save_plot('dCT_distribution_Q5');


figure;
plot(c_tip, P_total_Q5, '-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
yline(P_opt_quad, '--r', 'LineWidth', 2, 'Label', 'Desired Power', 'LabelHorizontalAlignment', 'left');
xlabel('Tip Chord ', 'FontSize', 12);
ylabel('Total Power [W]', 'FontSize', 12);
title('Total Power vs. Tip Chord Length', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
ylim([min(P_total_Q5), P_opt_quad + 20]);
save_plot('Power_tip_chord_Q5');


figure;
plot(c_tip, T_total_Q5, '-s', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
yline(T_opt_quad, '--g', 'LineWidth', 2, 'Label', 'Desired Thrust', 'LabelHorizontalAlignment', 'left');
xlabel('Tip Chord', 'FontSize', 12);
ylabel('Total Thrust [N]', 'FontSize', 12);
title('Total Thrust vs. Tip Chord Length', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
save_plot('Thrust_tip_chord_Q5');


% Find the index of c_tip_optimum
c_tip_optimum = 0.02;
idx_c_tip_opt = find(c_tip == c_tip_optimum);


% Retrieve values at c_tip_optimum
chord_E63_opt = chord_E63(:, idx_c_tip_opt);
radial_pos_E63_opt = radial_pos_E63(:, idx_c_tip_opt);
theta_E63_opt = theta_E63;
dP_E63_opt = dP_E63(:, idx_c_tip_opt);
dT_E63_opt = dT_E63(:, idx_c_tip_opt);
dCT_E63_opt = dCT_E63(:, idx_c_tip_opt);
T_total_Q5_opt = T_total_Q5(idx_c_tip_opt);
P_total_Q5_opt = P_total_Q5(idx_c_tip_opt);

% Find index corresponding to r/R = 0.75
idx_r_75 = find(y == 0.744848484848485);  % or == 0.75 if you know it's exact

% Extract the chord at r/R = 0.75 for the optimal c_tip
chord_E63_r_75 = 0.0455;



disp('Total thrust at c_tip_optimum:');
disp(T_total_Q5_opt);
disp('Total power at c_tip_optimum:');
disp(P_total_Q5_opt);

% Compute the percentage increase in thrust
percentage_increase_thrust = 100 * (T_total_Q5_opt / T_opt_quad - 1);

% Compute the percentage decrease in power
percentage_decrease_power = 100 * (1 - P_total_Q5_opt / P_opt_quad);


% Display the results
disp('Percentage Increase in Thrust from Thrust_opt_quad_final:');
disp(percentage_increase_thrust);
disp('Percentage Decrease in Power from Power_opt_quad_final:');
disp(percentage_decrease_power);

%% Task 6
m0 = 8.02; % Calculated from the slope of cl curve (near 0)
alpha_low = -2.2909; % in deg - AoA at lift = 0
density_blades = 74; % in kg/m^3
V_inf = 10;

CD_body = 0.4;

k_factor = T_opt_quad/(omega_opt_quad_q2^2);

beta = 0;       % Convert once

A_Q6 = pi*(R_opt_quad^2); % Area of each rotor


AR = linspace(6,12,4);
b_wing = linspace(0,4,100);
L_wing = zeros(length(b_wing), length(AR));
D_wing = zeros(length(b_wing), length(AR));
A_wing = zeros(length(b_wing), length(AR));
D_total_wing = zeros(length(b_wing), length(AR));
L_total_wing = zeros(length(b_wing), length(AR));
F_vertical_wing = zeros(length(b_wing), length(AR));
P_drag_wing = zeros(length(b_wing), length(AR));
P_total_Q6 = zeros(length(b_wing), length(AR));
m_total_Q6 = zeros(length(b_wing), length(AR));
m_blade = zeros(length(b_wing), length(AR));



for i = 1:length(b_wing)
    for k = 1:length(AR)
        chord_wing = b_wing(i) / AR(k);
        thickness_wing = 0.12*chord_wing;
        % Initialization for convergence loop
        diff_6 = 1;
        eps = 1e-6;
        
        % Initial guesses from previous optimal quad configuration
        m_old_6 = m_total_opt_quad;
        P_old_6 = P_opt_quad / 4;  % Per-rotor initial guess
        T_old_6 = T_opt_quad / 4;  % Per-rotor initial guess
        omega_q6 = omega_opt_quad_q2;
    
        while diff_6 > eps
            % Wing volume and blade mass for current span
            Volume_wing = b_wing(i) * chord_wing * thickness_wing;

            % Wing lift and drag only if wing exists
            if b_wing(i) == 0
                L_wing(i,k) = 0;
                D_wing(i,k) = 0;
                A_wing(i,k) = 0;
                D_total_wing(i,k) = 0;
                F_vertical_wing(i,k) = 0;
                m_blade(i,k) = 0;
            else
                [L_wing(i,k), D_wing(i,k)] = compute_lifting_line_wing(AR(k), b_wing(i), aoa_d_E63, m0, alpha_low, rho_mars, V_inf); 
                %[L_wing(i,k,n), D_wing(i,k,n)] =  Lifting_line_elliptical(AR(k), aoa_d_E63, alpha_low, target_cd, rho_mars, b_wing(i), chord_wing, V_inf);
   
                A_wing(i,k) = b_wing(i) * chord_wing;
                L_total_wing(i,k) = L_wing(i,k);
                D_profile = target_cd * (1/2) * rho_mars * (V_inf^2) *b_wing(i) * chord_wing;
                D_total_wing(i,k) = (D_wing(i,k)) + D_profile;
            end

            m_blade(i,k) = density_blades * Volume_wing;
            mass_battery_final = 0.5;
    
           
    
            m_total_Q6(i,k) = m_total_opt_quad + m_blade(i,k);
            
            % Induced velocity
            vh = sqrt(T_old_6 / (2 * rho_mars * A_Q6)); % Per rotor
            vi = induced_wind(vh, V_inf, beta);
    
            
            D_body_prop = 0.5*rho_mars* CD_body * (V_inf^2) *0.3;
            D_total = D_total_wing(i,k) + D_body_prop;
            % Power model
  
            [T_per_rotor_6, P_per_rotor_6, beta_new] = calculate_Power_Thrust_Q6(...
                (m_total_Q6(i,k)/4), grav_acc, rho_mars, R_opt_quad, ...
                chord_E63_r_75, 2, Cd_0, omega_q6, V_inf, vi, L_total_wing(i,k)/4, D_total/4);
            % Compute required total thrust to balance vertical + horizontal forces

            m_diff = abs(m_old_6 - m_total_Q6(i,k));
            P_diff = abs(P_old_6 - P_per_rotor_6);
            T_diff = abs(T_old_6 - T_per_rotor_6);
            beta_diff = abs(beta_new - beta);
            % Update convergence
            diff_6 = max([m_diff, P_diff, T_diff, beta_diff]);
            m_old_6 = m_total_Q6(i,k);
            P_old_6 = P_per_rotor_6;
            T_old_6 = T_per_rotor_6;
            beta = beta_new;
            omega_q6 = sqrt((4*T_per_rotor_6)/k_factor);
        end
    
        % Final outputs
  
        P_total_Q6(i,k) = P_per_rotor_6*4;
    end
end




figure;
hold on;

cmap = lines(length(AR));  % Use MATLAB's default line colors

for k = 1:length(AR)
    power_values = P_total_Q6(:, k);
    
    % Avoid division by zero
    if power_values(1) ~= 0
        normalized_power = power_values/power_values(1);
    else
        normalized_power = power_values;  % fallback if first value is zero
    end

    plot(b_wing, normalized_power, '-', 'LineWidth', 2, ...
         'Color', cmap(k,:), 'DisplayName', ['AR = ' num2str(AR(k))]);
end

xlabel('Wing Span b [m]');
ylabel('Normalized Power');
title('Effect of Wing Span and Aspect Ratio on Power');
legend('Location','best');
grid on;
save_plot('Power_with_wing_Q6');

figure;
hold on;

cmap = lines(length(AR));  % Use MATLAB's default line colors

for k = 1:length(AR)
    plot(b_wing, m_total_Q6(:,k), '-', 'LineWidth', 2, ...
         'Color', cmap(k,:), 'DisplayName', ['AR = ' num2str(AR(k))]);
end

xlabel('Wing Span b [m]');
ylabel('Mass [kg]');
title('Effect of Wing Span and Aspect Ratio on Mass');
legend('Location','best');
grid on;
save_plot('Mass_with_wing_Q6');

%% Task 6 - Part 2 - Corrected

density_blades_T2 = [10, 20, 30]; % in kg/m^3


beta_T2 = 0;  



AR_T2 = 150;
b_wing_T2 = linspace(0,40,100);
L_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
D_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
A_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
D_total_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
L_total_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
F_vertical_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
P_drag_wing_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
P_total_Q6_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
m_total_Q6_T2 = zeros(length(b_wing_T2), length(density_blades_T2));
m_blade_T2 = zeros(length(b_wing_T2), length(density_blades_T2));



for i = 1:length(b_wing_T2)
    for k = 1:length(density_blades_T2)
        chord_wing_T2 = b_wing_T2(i) / AR_T2;
        thickness_wing_T2 = 0.12*chord_wing_T2;
        % Initialization for convergence loop
        diff_6 = 1;
        eps = 1e-6;
        
        % Initial guesses from previous optimal quad configuration
        m_old_6_T2 = m_total_opt_quad;
        P_old_6_T2 = P_opt_quad / 4;  % Per-rotor initial guess
        T_old_6_T2 = T_opt_quad / 4;  % Per-rotor initial guess
        omega_q6_T2 = omega_opt_quad_q2;
    
        while diff_6 > eps
            % Wing volume and blade mass for current span
            Volume_wing_T2 = b_wing_T2(i) * chord_wing_T2 * thickness_wing_T2;

            % Wing lift and drag only if wing exists
            if b_wing_T2(i) == 0
                L_wing_T2(i,k) = 0;
                D_wing_T2(i,k) = 0;
                A_wing_T2(i,k) = 0;
                D_total_wing_T2(i,k) = 0;
                F_vertical_wing_T2(i,k) = 0;
                m_blade_T2(i,k) = 0;
            else
                [L_wing_T2(i,k), D_wing_T2(i,k)] = compute_lifting_line_wing(AR_T2, b_wing_T2(i), aoa_d_E63, m0, alpha_low, rho_mars, V_inf); 
                %[L_wing(i,k,n), D_wing(i,k,n)] =  Lifting_line_elliptical(AR(k), aoa_d_E63, alpha_low, target_cd, rho_mars, b_wing(i), chord_wing, V_inf);
   
                A_wing_T2(i,k) = b_wing_T2(i) * chord_wing_T2;
                L_total_wing_T2(i,k) = L_wing_T2(i,k);
                D_profile_T2 = target_cd * (1/2) * rho_mars * (V_inf^2) *b_wing_T2(i) * chord_wing_T2;
                D_total_wing_T2(i,k) = (D_wing_T2(i,k)) + D_profile_T2;
            end

            m_blade_T2(i,k) = density_blades_T2(k) * Volume_wing_T2;
            mass_battery_final = 0.5;
    
           
    
            m_total_Q6_T2(i,k) = m_total_opt_quad + m_blade_T2(i,k);
            
            % Induced velocity
            vh_T2 = sqrt(T_old_6_T2 / (2 * rho_mars * A_Q6));
            vi_T2 = induced_wind(vh_T2, V_inf, beta_T2);
    
            
            D_body_prop_T2 = 0.5*rho_mars* CD_body * (V_inf^2) *0.3;
            D_total_T2 = D_total_wing_T2(i,k) + D_body_prop_T2;
            % Power model
  
            [T_per_rotor_6_T2, P_per_rotor_6_T2, beta_new_T2] = calculate_Power_Thrust_Q6(...
                (m_total_Q6_T2(i,k)/4), grav_acc, rho_mars, R_opt_quad, ...
                chord_E63_r_75, 2, Cd_0, omega_q6_T2, V_inf, vi_T2, L_total_wing_T2(i,k)/4, D_total_T2/4);
            % Compute required total thrust to balance vertical + horizontal forces

            m_diff = abs(m_old_6_T2 - m_total_Q6_T2(i,k));
            P_diff = abs(P_old_6_T2 - P_per_rotor_6_T2);
            T_diff = abs(T_old_6_T2 - T_per_rotor_6_T2);
            beta_diff = abs(beta_new_T2 - beta_T2);
            % Update convergence
            diff_6 = max([m_diff, P_diff, T_diff, beta_diff]);
            m_old_6_T2 = m_total_Q6_T2(i,k);
            P_old_6_T2 = P_per_rotor_6_T2;
            T_old_6_T2 = T_per_rotor_6_T2;
            beta_T2 = beta_new_T2;
            omega_q6_T2 = sqrt((4*T_per_rotor_6_T2)/k_factor);
        end
    
        % Final outputs
  
        P_total_Q6_T2(i,k) = P_per_rotor_6_T2*4;
    end
end

figure;
hold on;
colors = lines(length(density_blades_T2));  % Use distinct colors for each density

for k = 1:length(density_blades_T2)
    plot(b_wing_T2, P_total_Q6_T2(:, k)/P_total_Q6_T2(1,k), 'LineWidth', 1.8, ...
        'Color', colors(k, :), ...
        'DisplayName', sprintf('Density = %.1f kg/m^3', density_blades_T2(k)));
end

% Add horizontal line at y = 0.9 without legend entry
h = yline(0.9, '--k', '10% reduction in Power', ...
    'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlabel('Wing Span b [m]');
ylabel('Normalized Power');
title('Total Power vs. Wing Span for Different Blade Densities');
legend('Location', 'best');
grid on;
save_plot('Power_with_wing_density_Q6(2)');



%% Task 6 - Part 3


density_blades_T3 = 74; % in kg/m^3
V_inf_T3 = linspace(0,12,26);

CD_body = 0.4;
eps_Q6 = 1e-3;
diff_Q6 = 1;
k_factor = T_opt_quad/(omega_opt_quad_q2^2);

beta = 0;       % Convert once




AR = 0;
b_wing_T3 = 0;
L_wing_T3 = zeros(length(V_inf_T3), 1);
D_wing_T3 = zeros(length(V_inf_T3), 1);
A_wing_T3 = zeros(length(V_inf_T3), 1);
D_total_wing_T3 = zeros(length(V_inf_T3), 1);
L_total_wing_T3 = zeros(length(V_inf_T3), 1);
F_vertical_wing_T3 = zeros(length(V_inf_T3), 1);
P_drag_wing_T3 = zeros(length(V_inf_T3), 1);
P_total_Q6_T3 = zeros(length(V_inf_T3), 1);
m_total_Q6_T3 = zeros(length(V_inf_T3), 1);
m_blade_T3 = zeros(length(V_inf_T3), 1);
beta_new_T3 = zeros(length(V_inf_T3), 1);



for i = 1:length(V_inf_T3)

        % Initialization for convergence loop
        diff_6 = 1;
        eps = 1e-6;
        
        % Initial guesses from previous optimal quad configuration
        m_old_6_T3 = m_total_opt_quad;
        P_old_6_T3 = P_opt_quad / 4;  % Per-rotor initial guess
        T_old_6_T3 = T_opt_quad / 4;  % Per-rotor initial guess
        omega_q6_T3 = omega_opt_quad_q2;
    
        while diff_6 > eps


            % Wing lift and drag only if wing exists
            if b_wing_T3 == 0
                L_wing_T3(i) = 0;
                D_wing_T3(i) = 0;
                A_wing_T3(i) = 0;
                D_total_wing_T3(i) = 0;
                F_vertical_wing_T3(i) = 0;
                m_blade_T3(i) = 0;
            else
                [L_wing_T3(i), D_wing_T3(i)] = compute_lifting_line_wing(AR, b_wing_T3, aoa_d_E63, m0, alpha_low, rho_mars, V_inf_T3(i)); 
                %[L_wing(i,k,n), D_wing(i,k,n)] =  Lifting_line_elliptical(AR(k), aoa_d_E63, alpha_low, target_cd, rho_mars, b_wing(i), chord_wing, V_inf);
   
                A_wing_T3(i) = b_wing_T3 * chord_wing;
                L_total_wing_T3(i) = L_wing_T3(i);
                D_profile = target_cd * (1/2) * rho_mars * (V_inf_T3(i)^2) *b_wing_T3 * chord_wing;
                D_total_wing_T3(i) = (D_wing_T3(i)) + D_profile;
            end


            mass_battery_final = 0.5;
    
           
    
            m_total_Q6_T3(i) = m_total_opt_quad + m_blade_T3(i);
            
            % Induced velocity
            vh_T3 = sqrt(T_old_6_T3 / (2 * rho_mars * A_Q6));
            vi_T3 = induced_wind(vh_T3, V_inf_T3(i), beta);
    
            
            D_body_prop_T3 = 0.5*rho_mars* CD_body * (V_inf_T3(i)^2) *0.3;
            D_total_T3 = D_total_wing_T3(i) + D_body_prop_T3;
            % Power model
  
            [T_per_rotor_6_T3, P_per_rotor_6_T3, beta_new_T3(i)] = calculate_Power_Thrust_Q6(...
                (m_total_Q6_T3(i)/4), grav_acc, rho_mars, R_opt_quad, ...
                chord_E63_r_75, 2, Cd_0, omega_q6_T3, V_inf_T3(i), vi_T3, L_total_wing_T3(i)/4, D_total_T3/4);
            % Compute required total thrust to balance vertical + horizontal forces

            m_diff = abs(m_old_6_T3 - m_total_Q6_T3(i));
            P_diff = abs(P_old_6_T3 - P_per_rotor_6_T3);
            T_diff = abs(T_old_6_T3 - T_per_rotor_6_T3);
            beta_diff = abs(beta_new_T3(i) - beta);
            % Update convergence
            diff_6 = max([m_diff, P_diff, T_diff, beta_diff]);
            m_old_6_T3 = m_total_Q6_T3(i);
            P_old_6_T3 = P_per_rotor_6_T3;
            T_old_6_T3 = T_per_rotor_6_T3;
            beta = beta_new_T3(i);
            omega_q6_T3 = sqrt((4*T_per_rotor_6_T3)/k_factor);
        end
    
        % Final outputs
  
        P_total_Q6_T3(i) = P_per_rotor_6_T3*4;

end

% Ensure your vectors/matrices are defined properly
% V_inf_T3 should be a vector of free-stream velocities
% P_total_Q6_T3 and beta_new_T3 should be vectors or have size matching V_inf_T3

figure;
yyaxis left
plot(V_inf_T3, P_total_Q6_T3, '-o', 'LineWidth', 2);
ylabel('Total Power P_{total} (W)');
xlabel('Free-stream Velocity V_{\infty} (m/s)');
grid on;
title('Power Consumption and \beta vs V_{\infty}');


yyaxis right
plot(V_inf_T3, rad2deg(beta_new_T3), '--s', 'LineWidth', 2);
ylabel('\beta Angle (deg)');
ylim([0, 1.2*max(rad2deg(beta_new_T3))]);

legend('Total Power (W)', '\beta Angle (deg)', 'Location', 'best');
save_plot('Power_beta_final_Q6(3)');

%% Functions

function [m_propeller, m_control, m_computer, m_no_fuselage, m_fuselage, m_total] = calc_mass(n_blades, R_new, P_new, P_old, n_rotors, new_batt_mass, m_payload)
    R_old = 0.605;  % m

    m_fuz_ing = 0.3;  % kg
    m_no_fuz_ing = 1.5;  % kg
    m_ing = m_no_fuz_ing + m_fuz_ing;  % kg

    m_propeller = (0.07/4 * n_blades * R_new/R_old)*n_rotors;  % kg
    m_control = (0.25/2 * P_new / P_old)*n_rotors;  % kg
    m_computer = 1;  % kg

    m_no_fuselage = m_computer + m_control + m_propeller + new_batt_mass ;  % kg
    m_fuselage = ((m_no_fuselage)) * (m_fuz_ing / m_no_fuz_ing);  % kg
    m_total = m_no_fuselage + m_fuselage + m_payload;


end


function [T, P_tot] = calculate_Power_Thrust(m_total, grav_acc, nr_prop, rho, R, needed_c_div_R, n_blades, Cd_0, omega)
    % Compute thrust required for half the mass (assuming 2 rotors)
    gamma = 1.15;
    T = m_total * grav_acc / nr_prop;
    
    % Ideal power for one rotor
    P_ideal = sqrt(T^3) / sqrt(2 * rho * pi * R^2);
    finally_c = needed_c_div_R * (R);
    
    % Profile power losses (induced + profile drag)
    P_0 = (1/8) * rho * finally_c * n_blades * Cd_0 * omega^3 * R^4;
    
    % Total power for one rotor with losses
    P_tot = P_ideal * gamma + P_0;
end

function [Thrust_per_rotor, P_tot, beta] = calculate_Power_Thrust_Q6(mass, grav_acc, rho, R, chord, n_blades, Cd_0, omega, V_inf,  vi, Lift_wing, Drag_total)
    % Compute thrust required for half the mass (assuming 2 rotors)
    gamma = 1.15;
    Weight_per_rotor = mass*grav_acc;
    % Lift_req = Weight_per_rotor - Lift_wing;
    %beta = atan(Drag_total / (Weight_per_rotor));
    beta = atan(Drag_total / (Weight_per_rotor - Lift_wing));
    %Thrust_per_rotor = Lift_req/cos(beta);
    Thrust_per_rotor = sqrt(((Weight_per_rotor - Lift_wing)^2) + (Drag_total^2));
    
   
    

    % Ideal power for one rotor
    P_ideal = Thrust_per_rotor*((V_inf*sin(beta)) + vi);
    finally_c = chord;
    
    % Profile power losses (induced + profile drag)
    P_0 = (1/8) * rho * finally_c * n_blades * Cd_0 * omega^3 * R^4;
    
    % Total power for one rotor with losses
    P_tot = P_ideal * gamma + P_0;
end

function chord = chord_airfoil(c_tip, r, R)
    chord = c_tip .* (1 + 5 * (1 - r./R));  % root is ~6× tip
end








function vi = induced_wind(v_H, V_inf, beta)
   eqn = @(vi) vi - ((v_H^2) / (sqrt(((V_inf*cos(beta))^2) + ((V_inf*sin(beta) + vi)^2))));

    % Solve using fsolve
    vi0 = 1;  % Initial guess
    options = optimoptions('fsolve','Display','off');
    vi = fsolve(eqn, vi0, options);
end







function [L_wing, Di_wing] = compute_lifting_line_wing(AR,b_span, alpha_d, m0, alpha_low, rho_mars, V_inf) 

    AoA_Q2 = alpha_d;
    n = 50;
    C = b_span/AR;
    theta = linspace(0.01, pi, 180);
    A_interm = zeros(n,1);
    A_interm(:) = compute_An_q3(n, b_span, m0, C, AoA_Q2, alpha_low, theta);

    Gamma_inter_q3 = zeros(n, length(theta));
    Gamma_q3 = zeros(length(theta),1);
    for i = 1:length(theta)  % Loop over rows of theta_q3
        Gamma_sum_q3 = 0;
        for j = 1:n
            Gamma_inter_q3(j, i) = A_interm(j) * sin(j * theta(i));  % Fixed indexing
            Gamma_sum_q3 = Gamma_inter_q3(j, i) + Gamma_sum_q3;
        end
            
        Gamma_q3(i) = Gamma_sum_q3 * 2 * b_span *V_inf;  % Corrected indexing
    end
    
    CL = pi*AR*A_interm(1);
    value = zeros(length(n),1);
    for i = 2:n
        value(i) = i * ((A_interm(i) / A_interm(1))^2);
    end
    delta = sum(value(2:n));  % sum from i = 2 to n
    CDi = ((CL^2)/(pi*AR))*(1+delta);
    L_wing = CL*(0.5)*rho_mars*(b_span*C)*(V_inf^2);
    Di_wing = CDi*(0.5)*rho_mars*(b_span*C)*(V_inf^2);
   
end

function An = compute_An_q3(n, b, m0, C, alpha, alpha_L0, theta)
     % SOLVE_COEFFICIENTS Solves for the first n coefficients A_n
    % Inputs:
    %   n        - Number of coefficients to compute
    %   b        - Given constant parameter
    %   m0       - Function handle or vector of m0(theta)
    %   C        - Function handle or vector of C(theta)
    %   alpha    - Function handle or vector of alpha(theta)
    %   alpha_L0 - Function handle or vector of alpha_L0(theta)
    %   theta    - Vector of theta values for numerical evaluation

    
    % Construct the system of equations
    S = zeros(length(theta), n);
    T = zeros(length(theta), n);
    
    for k = 1:n
        for j = 1:length(theta)
            if theta(j) == 0
                theta(j) = 1e-6;
            end
            % if theta(j) == pi
            %     theta(j) = 0.01;
            % end
            S(j, k) = (-4*b*sin(k*theta(j))) / (m0 * C);
            T(j, k) = (k*sin(k*theta(j))/sin(theta(j)));
        end
    end
    
    % Right-hand side of the equation
    RHS = zeros(length(theta),1);
    for i = 1:length(theta)
        RHS(i) = -deg2rad(alpha - alpha_L0);
    end

    
    % Solve the system
    A_matrix = S - T;
    An = A_matrix \  RHS;

end


function [dT_BE, dP, dCT, aoa, cl, cd]  = calculate_BEM_final(omega, r, theta, aoa_E63, cl_E63, cd_E63, Num_blades, rho, chord, R)
    diff = 1;
    eps = 1e-4;
    vi = 0.01;
    max_iter = 1000;
    iter = 0;

    while diff > eps && iter < max_iter
        V_rel = sqrt((vi^2) + ((omega * r)^2));
        phi = acos((omega * r)/(V_rel));
        aoa = deg2rad(theta) - phi;
        aoa_deg = theta - rad2deg(phi);  % Get AOA in degrees
        cl = interp1(aoa_E63, cl_E63, aoa_deg, 'linear', 'extrap');
        cd = interp1(aoa_E63, cd_E63, aoa_deg, 'linear', 'extrap');


        f_inter = (Num_blades/2)*((R-r)/(r*sin(phi)));
        F = (2/pi)*acos(exp(-f_inter));
        dT_BE = (1/2)*Num_blades*rho*chord*(V_rel^2)*((cl*cos(phi)) - (cd*sin(phi)));
        dT_MOM = 4*pi*rho*(vi^2)*r;
        dP = (1/2)*Num_blades*omega*rho*chord*(V_rel^2)*((cl*sin(phi)) + (cd*cos(phi)))*r;
        V_tip = omega*R;
        A = pi*(R^2);
        dCT = dT_BE / (rho*A*(V_tip^2));
        vi_new = vi + 0.5*(dT_BE - (F*dT_MOM));  % enforce positive vi
        diff = abs(vi_new - vi);
        if diff > eps
            vi = vi_new;
            iter = iter +1;
        end
    end
end

    

function theta = twist_distribution(r, R, theta_0, theta_tip, n)
    theta = theta_0 * (1 - (r./R).^n) + theta_tip;
end

function save_plot(file_name)
    % Ajusta el tamaño para documentos LaTeX (~6x4 pulgadas)
    set(gcf, 'Units', 'inches', 'Position', [1 1 6 4]);
    set(gca, 'FontSize', 11); % tamaño consistente para LaTeX
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 4]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 6 4]);

    % Guarda como PDF vectorial de alta calidad
    print(gcf, file_name, '-dpdf');
end