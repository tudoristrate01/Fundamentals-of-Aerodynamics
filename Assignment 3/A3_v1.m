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
old_batt_capaicty = 10;  % Wh
new_batt_capaicty = 20;  % Wh

n_blades = [2, 3, 4];

m_no_fuselage_interm = zeros(length(R_new), 2, 3);
m_fuselage_interm = zeros(length(R_new), 2, 3);

m_no_fuselage = zeros(length(R_new), 2, 3);
m_fuselage = zeros(length(R_new), 2, 3);

P_new = zeros(length(R_new), 2, 3);

for k  = 1 : length(n_blades)
    for j = 1 : length(nr_prop)
        for i = 1 : length(R_new)
            [~, ~, ~, m_no_fuselage_interm(i, j, k), m_fuselage_interm(i, j, k)] = calc_mass(n_blades(k), R_new(i), P_tot, P_tot, nr_prop(j));
    
            T_2 = m_fuselage_interm(i, j, k) * grav_acc / 2;
            P_ideal_2 = sqrt(T_2^3) / sqrt(2*rho*pi * R_new(i)^2);  % for 1 rotor
    
            P_0_2 = 1/8 * rho * finally_c * n_blades(k) * Cd_0 * omega^3 * R_new(i)^4;
    
            P_tot_2 = P_ideal_2 * gamma + P_0_2;
            omega_new = (P_tot_2/T_2)/R_new(i);
    
            [~, ~, ~, m_no_fuselage(i, j, k), m_fuselage(i, j, k)] = calc_mass(n_blades(k), R_new(i), P_tot_2, P_tot, nr_prop(j));
            T_new = m_fuselage(i, j, k) * grav_acc / 2;
            P_ideal_new = sqrt(T_new^3) / sqrt(2*rho*pi * R_new(i)^2);  % for 1 rotor
    
            P_0_new = 1/8 * rho * finally_c * n_blades(k) * Cd_0 * omega_new^3 * R_new(i)^4;
    
            P_new(i, j, k) = P_ideal_new * gamma + P_0_new;
        end
    end
end

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


%% Question 3

% For bicopter (nr_prop(1))
[P_bi_min, idx_bi] = min(P_new(:, 1, 1));  % j = 1, k = 1
R_opt_bi = R_new(idx_bi);
P_opt_bi = P_bi_min;

% For quadcopter (nr_prop(2))
[P_quad_min, idx_quad] = min(P_new(:, 2, 1));  % j = 2, k = 1
R_opt_quad = R_new(idx_quad);
P_opt_quad = P_quad_min;



%% Functions

function [m_propeller, m_control, m_computer, m_no_fuselage, m_fuselage] = calc_mass(n_blades, R_new, P_new, P_old, n_rotors)
    R_old = 0.605;  % m

    m_fuz_ing = 0.3;  % kg
    m_no_fuz_ing = 1.5;  % kg
    m_ing = m_no_fuz_ing + m_fuz_ing;  % kg

    m_propeller = 0.07/4 * n_blades * R_new/R_old;  % kg
    m_control = 0.25/n_rotors * P_new / P_old;  % kg
    m_computer = 1;  % kg

    m_no_fuselage = m_computer + m_control + m_propeller;  % kg
    m_fuselage = m_no_fuselage * (m_fuz_ing / m_no_fuz_ing);  % kg

end