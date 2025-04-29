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

P_tot = gamma * P_ideal + P_0;


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

omega_2 = zeros(length(R_new), 2, 3);

P_new = zeros(length(R_new), 2, 3);

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
                [~, ~, ~, m_no_fuselage_interm(i, j, k), m_fuselage_interm(i, j, k), m_total(i,j,k)] = ...
                    calc_mass(n_blades(k), R_new(i), P_old, P_tot, nr_prop(j), new_batt_mass);
                [P_per_rotor] = calculate_Power_Thrust(m_total(i,j,k), grav_acc, nr_prop(j), rho, R_new(i), needed_c_div_R, n_blades(k), Cd_0, omega_old);
                P_new(i,j,k) = P_per_rotor * nr_prop(j);
                omega_2(i,j,k) = tip_speed_ing * R_new(i);
                diff = abs(m_old - m_total(i,j,k));
                if diff > eps
                    m_old = m_total(i,j,k);
                    omega_old = omega_2(i,j,k);
                    P_old = P_new(i,j,k);
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
omega_opt_quad = omega_2(idx_bi, 2, 1);
m_total_opt_quad = m_total(idx_quad, 2, 1);
flight_time_quad = (new_batt_capaicty/(P_opt_quad))/60; % In minutes

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
fprintf('Optimal Angular Velocity (omega_opt_quad): %.3f rad/s\n', omega_opt_quad);
fprintf('Total Mass (m_total_opt_quad): %.3f kg\n', m_total_opt_quad);
fprintf('Estimated Flight Time: %.2f minutes\n', flight_time_quad);


% Plot for Bicopter
figure(1);
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

% Plot for Quadcopter
figure(2);
hold on;
for i = 1 : 3
    plot(R_new, P_new(:, 2, i), 'LineWidth', 1.5);
end

% Mark the optimal point for quadcopter (2 blades, j=2, k=1)
plot(R_opt_quad, P_opt_quad, 'go', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Radius (m)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
grid on
title('Quadcopter')
legend({'2 blades', '3 blades', '4 blades', 'Optimal (2 blades)'}, 'Location', 'best');
hold off;

%% Task 3

each_battery_mass = 0.047; % kg
each_battery_capacity = (10/6)*3600; % J
n_add_battery = linspace(0,300,301);
n_battery = 6 + n_add_battery;

m_no_fuselage_3 = zeros(length(n_battery),1);
m_fuselage_3 = zeros(length(n_battery),1);
m_total_3 = zeros(length(n_battery),1);
Power_3 = zeros(length(n_battery),1);
flight_time_3 = zeros(length(n_battery),1);


for i = 1:length(n_battery)
    diff_3 = 20;
    eps = 1e-6;
    m_old_3 = m_total_opt_quad;
    P_old_3 = P_opt_quad;
    while diff_3 > eps
        m_battery = n_battery(i)*each_battery_mass;
        cap_battery = n_battery(i)*each_battery_capacity;
        [~, ~, ~, m_no_fuselage_3, m_fuselage_3, m_total_3(i)] = calc_mass(2, R_opt_quad, P_old_3, P_opt_quad, 2, m_battery);
        [Power_3(i)] = calculate_Power_Thrust(m_total_3(i), grav_acc, 4, rho, R_opt_quad, needed_c_div_R, 2, Cd_0, omega_opt_quad);
        diff_3 = abs(m_old_3 - m_total_3(i));
        if diff_3 > eps
            m_old_3 = m_total_3(i);
            P_old_3 = Power_3(i);
        end


        flight_time_3(i) = cap_battery/(Power_3(i)*4);
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


% Display the result
fprintf('\n========= Quadcopter Flight Time for Optimal no. of %d Batteries =========\n', max_pow_batt);
fprintf('Flight Time: %.2f minutes\n', max_flight_time_quad / 60);
fprintf('Total Batteries Used: %.0f\n', max_pow_batt);

% Display the result
fprintf('\n========= Quadcopter Flight Time for Limit of %d Batteries =========\n', limit_batt);
fprintf('Flight Time: %.2f minutes\n', flight_time_limit_batt / 60);
fprintf('Total Batteries Used: %.0f\n', limit_batt);


figure(3)
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



%% Functions

function [m_propeller, m_control, m_computer, m_no_fuselage, m_fuselage, m_total] = calc_mass(n_blades, R_new, P_new, P_old, n_rotors, new_batt_mass)
    R_old = 0.605;  % m

    m_fuz_ing = 0.3;  % kg
    m_no_fuz_ing = 1.5;  % kg
    m_ing = m_no_fuz_ing + m_fuz_ing;  % kg

    m_propeller = 0.07/4 * n_blades * R_new/R_old;  % kg
    m_control = 0.25/n_rotors * P_new / P_old;  % kg
    m_computer = 1;  % kg

    m_no_fuselage = m_computer + m_control + m_propeller + new_batt_mass +2 ;  % kg
    m_fuselage = ((m_no_fuselage*n_rotors)) * (m_fuz_ing / m_no_fuz_ing);  % kg
    m_total = m_no_fuselage + m_fuselage;


end


function [P_tot] = calculate_Power_Thrust(m_total, grav_acc, nr_prop, rho, R, needed_c_div_R, n_blades, Cd_0, omega)
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
