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
R_new = linspace(0.3, 0.7, 41); 
new_batt_mass = 0.5;  % kg
old_batt_capaicty = 10*3600;  % J
new_batt_capaicty = 20*3600;  % J

n_blades = [2, 3, 4];

m_no_fuselage_interm = zeros(length(R_new), 2, 3);
m_fuselage_interm = zeros(length(R_new), 2, 3);

m_no_fuselage = zeros(length(R_new), 2, 3);
m_fuselage = zeros(length(R_new), 2, 3);

omega_2 = zeros(length(R_new), 2, 3);

P_new = zeros(length(R_new), 2, 3);

for k = 1:length(n_blades)
    for j = 1:length(nr_prop)
        for i = 1:length(R_new)

            % Initialize diff to enter the while loop
            diff = inf;

            while diff > 1e-1
                % First mass estimation with given total power (P_tot)
                [~, ~, ~, m_no_fuselage_interm(i, j, k), m_fuselage_interm(i, j, k)] = ...
                    calc_mass(n_blades(k), R_new(i), P_tot, P_tot, nr_prop(j), new_batt_mass);
                
                % Compute thrust required for half the mass (assuming 2 rotors)
                T_2 = m_fuselage_interm(i, j, k) * grav_acc / 2;

                % Ideal power for one rotor
                P_ideal_2 = sqrt(T_2^3) / sqrt(2 * rho * pi * R_new(i)^2);

                % Profile power losses (induced + profile drag)
                P_0_2 = (1/8) * rho * finally_c * n_blades(k) * Cd_0 * omega^3 * R_new(i)^4;

                % Total power for one rotor with losses
                P_tot_2 = P_ideal_2 * gamma + P_0_2;

                % New angular velocity estimate
                omega_2(i, j, k) = tip_speed_ing * R_new(i);
                omega_new = omega_2(i, j, k);

                % Updated mass estimation with refined power value
                [~, ~, ~, m_no_fuselage(i, j, k), m_fuselage(i, j, k)] = ...
                    calc_mass(n_blades(k), R_new(i), P_tot_2, P_tot, nr_prop(j), new_batt_mass);
                
                % Recalculate thrust with updated mass
                T_new = m_fuselage(i, j, k) * grav_acc / 2;

                % Recalculate ideal power for one rotor
                P_ideal_new = sqrt(T_new^3) / sqrt(2 * rho * pi * R_new(i)^2);

                % Recalculate profile power with new omega
                P_0_new = (1/8) * rho * finally_c * n_blades(k) * Cd_0 * omega_new^3 * R_new(i)^4;

                % Recalculate total power
                P_new(i, j, k) = P_ideal_new * gamma + P_0_new;

                % Update convergence criteria
                diff = abs(m_fuselage_interm(i, j, k) - m_fuselage(i, j, k));
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
m_fuselage_opt_bi = m_fuselage(idx_bi, 1, 1);  % j = 1, k = 1
flight_time_bi = new_batt_capaicty/P_opt_bi;

% For quadcopter (nr_prop(2)) - 2 Blades shows optimal performance
[P_quad_min, idx_quad] = min(P_new(:, 2, 1));  % j = 2, k = 1
R_opt_quad = R_new(idx_quad);
P_opt_quad = P_quad_min;
omega_opt_quad = omega_2(idx_bi, 2, 1);
m_fuselage_opt_quad = m_fuselage(idx_quad, 2, 1);
flight_time_quad = new_batt_capaicty/P_opt_quad;


% Plots
figure(1);
hold on;
for i = 1 : 3
    plot(R_new, P_new(:, 1, i), 'LineWidth', 1.5);
end
xlabel('Radius', 'FontSize', 12);
ylabel('Power', 'FontSize', 12);
grid on
title('Bicopter')
legend({'2 blades', '3 blades', '4 blades'}, 'Location', 'best');
hold off;

figure(2);
hold on;
for i = 1 : 3
    plot(R_new, P_new(:, 2, i), 'LineWidth', 1.5);
end
xlabel('Radius', 'FontSize', 12);
ylabel('Power', 'FontSize', 12);
grid on
title('Quadcopter')
legend({'2 blades', '3 blades', '4 blades'}, 'Location', 'best');
hold off;

%% Task 3

each_battery_mass = 0.047; % kg
each_battery_capacity = (10/6)*3600; % J
n_add_battery = linspace(1,60,61);
n_battery = 6 + n_add_battery;

omega_3 = zeros(length(n_battery),1);
m_no_fuselage_3 = zeros(length(n_battery),1);
m_fuselage_3 = zeros(length(n_battery),1);
Power_3 = zeros(length(n_battery),1);
flight_time_3 = zeros(length(n_battery),1);


for i = 1:length(n_battery)
    m_battery = n_battery(i)*each_battery_mass;
    cap_battery = n_battery(i)*each_battery_capacity;
    [~, ~, ~, m_no_fuselage_interm, m_fuselage_interm] = calc_mass(2, R_opt_quad, P_opt_quad, P_tot, 2, m_battery);
    
    T_3 = m_fuselage_interm* grav_acc / 2;
    P_ideal_3 = sqrt(T_3^3) / sqrt(2*rho*pi * R_opt_quad^2);  % for 1 rotor

    P_0_3 = 1/8 * rho * finally_c * 2 * Cd_0 * omega_opt_quad^3 * R_opt_quad^4;

    P_tot_3 = P_ideal_3 * gamma + P_0_3;
    omega_3(i) = (P_tot_3/T_3)/R_opt_quad;
    omega_new_3 = omega_3(i);


    [~, ~, ~, m_no_fuselage_3(i), m_fuselage_3(i)] = calc_mass(2, R_opt_quad, P_tot_3, P_tot, 2, m_battery);
    T_new_3 = m_fuselage_3(i) * grav_acc / 2;
    P_ideal_new_3 = sqrt(T_new_3^3) / sqrt(2*rho*pi * R_opt_quad^2);  % for 1 rotor

    P_0_new_3 = 1/8 * rho * finally_c * 2 * Cd_0 * omega_new_3^3 * R_opt_quad^4;

    Power_3(i) = P_ideal_new_3 * gamma + P_0_new_3;
    flight_time_3(i) = cap_battery/Power_3(i);
end

% Find the highest flight time and its corresponding n_battery
[max_flight_time_quad, idx] = max(flight_time_3);  % idx gives the index of the max flight time
max_add_battery_num_quad = n_add_battery(idx);  % Get the corresponding n_battery

% Display the result
disp(['Highest flight time: ', num2str(max_flight_time_quad/60), ' minutes']);
disp(['Corresponding n_battery: ', num2str(max_add_battery_num_quad)]);

figure(3)
plot(n_add_battery, flight_time_3/60);
xlabel('No. of batteries');
ylabel('Flight time (min)');
xline(42, 'r', 'LineWidth', 2);  % Adds a red vertical line at x = 42
grid on;
legend( 'Flight time - quadcopter', 'Maximum number of batteries - 42', 'Location', 'southwest');
title('Flight time vs. No. of batteris - Quadcopter');


%% Functions

function [m_propeller, m_control, m_computer, m_no_fuselage, m_fuselage] = calc_mass(n_blades, R_new, P_new, P_old, n_rotors, new_batt_mass)
    R_old = 0.605;  % m

    m_fuz_ing = 0.3;  % kg
    m_no_fuz_ing = 1.5;  % kg
    m_ing = m_no_fuz_ing + m_fuz_ing;  % kg

    m_propeller = 0.07/4 * n_blades * R_new/R_old;  % kg
    m_control = 0.25/n_rotors * P_new / P_old;  % kg
    m_computer = 1;  % kg

    m_no_fuselage = m_computer + m_control + m_propeller + new_batt_mass ;  % kg
    m_fuselage = m_no_fuselage * (m_fuz_ing / m_no_fuz_ing);  % kg

end


