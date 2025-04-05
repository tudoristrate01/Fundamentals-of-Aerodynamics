clc
clear
close all

%% Question 1 Calculations

Re = 6e6;
max_thicc = 0.15;
max_thicc_from_LE = 30.9;  % where the max thickenss is relative to the the leading edge
max_camber = 0.04;
max_camber_from_LE = 0.402;

fileID = fopen('cd.txt', 'r');

for i = 1:12
    fgetl(fileID);
end

data = textscan(fileID, '%f %f %f %f %f %f %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fileID);

CD_f = data{3};

alpha_low = -4.2; % Alpha for C_l closest to 0

Aspect_Ratio = [4, 6, 8, 10, inf];
AoA = linspace(-6, 10, 17);

C_l = zeros(length(Aspect_Ratio), length(AoA));
Cd_i = zeros(length(Aspect_Ratio), length(AoA));
CD = zeros(length(Aspect_Ratio), length(AoA));
alpha_i_rad = zeros(length(Aspect_Ratio), length(AoA));

for i = 1:length(Aspect_Ratio)
    for j = 1:length(AoA)
        [C_l(i, j), Cd_i(i, j), alpha_i_rad(i, j)] = calculate_Q1(Aspect_Ratio(i), AoA(j), alpha_low);
        CD(i, j) = Cd_i(i, j) + CD_f(j);
    end
end

alpha_i = rad2deg(alpha_i_rad);

%% Figures Q1

figure (1)
hold on
for i = 1:length(Aspect_Ratio)
    plot (AoA, C_l(i, :), LineWidth=3)
end
hold off
grid on
legend('ar4', 'ar6', 'ar8', 'ar10', 'ar inf')
xlabel('Angle of attack [Deg]')
ylabel('Lift coefficient [-]')
title('Lift coefficient vs Angle of attack')

figure (2)
hold on
for i = 1:length(Aspect_Ratio)
    plot (AoA, Cd_i(i, :), LineWidth=3)
end
hold off
grid on
legend('ar4', 'ar6', 'ar8', 'ar10', 'ar inf')
xlabel('Angle of attack [Deg]')
ylabel('Induced Drag coefficent [-]')
title('Induced Drag coefficent vs Angle of attack')

figure (3)
hold on
for i = 1:length(Aspect_Ratio)
    plot (AoA, CD(i, :), LineWidth=3)
end
hold off
grid on
legend('ar4', 'ar6', 'ar8', 'ar10', 'ar inf')
xlabel('Angle of attack [Deg]')
ylabel('Drag coefficent [-]')
title('Drag coefficent vs Angle of attack')

figure (4)
hold on
for i = 1:length(Aspect_Ratio)
    plot (AoA, alpha_i(i, :), LineWidth=3)
end
hold off
grid on
legend('ar4', 'ar6', 'ar8', 'ar10', 'ar inf')
xlabel('Angle of attack [Deg]')
ylabel('Induced angle of attack [Deg]')
title('Induced angle of attack vs Angle of attack')

close all   %  DELETE THIS ROW ===================================================================

%% Question 2 Calculations

AoA_Q2 = [0, 5, 10];

m0 = 6.5958;
C = 1;

n = 50;
B = Aspect_Ratio * C;
theta = linspace(0.1, pi, 180);
A_interm = zeros(n, length(AoA_Q2), length(B));

for i = 1:length(AoA_Q2)
    for j = 1:length(B)
        A_interm(:, i, j) = compute_An(n, B(j), m0, C, AoA_Q2(i), alpha_low, theta);
    end
end

An_plot = zeros(length(theta),1);

for i = 1:length(theta)
    An_plot(i) = A_interm(5,1,1)*sin(5*theta(i));
end

figure (5)
plot(theta/pi, An_plot);

%  Local properties from Slide 26

alpha_i_theta = zeros(length(theta), length(AoA_Q2), length(B));

for i = 1:length(AoA_Q2)
    for j = 1:length(B)
        An = A_interm(:, i, j);
        sum_term = zeros(size(theta));
        
        for n = 1:n
            sum_term = sum_term + n * An(n) * sin(n * theta) ./ sin(theta);
        end
        
        alpha_i_theta(:, i, j) = sum_term;
    end
end

Cl = zeros(length(theta), length(AoA_Q2), length(B));

for i = 1:length(AoA_Q2)
    for j = 1:length(B)
        An = A_interm(:, i, j);
        sum_term = zeros(size(theta));
        
        for n = 1:n
            sum_term = sum_term + An(n) * sin(n * theta);
        end
        
        Cl(:, i, j) = (4 * B(j) ./ C) .* sum_term;
    end
end

Cdi = zeros(length(theta), length(AoA_Q2), length(B));

for i = 1:length(AoA_Q2)
    for j = 1:length(B)
        An = A_interm(:, i, j);
        sum1 = zeros(size(theta));
        sum2 = zeros(size(theta));
        
        for n = 1:n
            sum1 = sum1 + An(n) * sin(n * theta);
        end
        
        for k = 1:n
            sum2 = sum2 + k * An(k) * (sin(k * theta) ./ sin(theta));
        end
        
        Cdi(:, i, j) = (4 * B(j) ./ C) .* sum1 .* sum2;
    end
end

% Convert theta to spanwise position y/b
y_b = cos(theta); % spanwise location normalized by half-span

% Define colors for better visualization
colors = lines(length(Aspect_Ratio));

% Plot results for each AoA
for i = 1:length(AoA_Q2)
    figure;
    hold on;
    grid on;
    
    for j = 1:length(Aspect_Ratio)
        plot(y_b, rad2deg(alpha_i_theta(:, i, j)), 'Color', colors(j, :), 'LineWidth', 2);
    end
    
    xlabel('Spanwise Position (y/b)');
    ylabel('\alpha_i (\theta) [deg]');
    title(['Induced Angle of Attack Distribution at AoA = ', num2str(AoA_Q2(i)), '°']);
    legend(arrayfun(@(x) sprintf('AR = %g', x), Aspect_Ratio, 'UniformOutput', false), 'Location', 'Best');
    
    hold off;
end

% Define line styles for better visualization
lineStyles = {'-', '--', ':', '-.', '-'};

% Plot Cl vs AoA
figure;
hold on;
for j = 1:length(Aspect_Ratio)
    plot(AoA_Q2, squeeze(mean(Cl(:, :, j), 1)), lineStyles{j}, 'LineWidth', 1.5, 'DisplayName', ['AR = ', num2str(Aspect_Ratio(j))]);
end
hold off;
grid on;
xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_l)');
title('Lift Coefficient vs. Angle of Attack');
legend('Location', 'best');
set(gca, 'FontSize', 12);

% Plot Cd_i vs AoA
figure;
hold on;
for j = 1:length(Aspect_Ratio)
    plot(AoA_Q2, squeeze(mean(Cdi(:, :, j), 1)), lineStyles{j}, 'LineWidth', 1.5, 'DisplayName', ['AR = ', num2str(Aspect_Ratio(j))]);
end
hold off;
grid on;
xlabel('Angle of Attack (degrees)');
ylabel('Induced Drag Coefficient (C_{d_i})');
title('Induced Drag Coefficient vs. Angle of Attack');
legend('Location', 'best');
set(gca, 'FontSize', 12);


%% Question 3

AR_3 = 6;
AoA_3 = deg2rad(6);
Taper_ratios = [0.2, 0.4, 0.6, 0.8, 1];



%% Functions

disp('Code finished!')

function [C_l, Cd_i, alpha_i_rad] = calculate_Q1(AR, alpha, alpha_low)  % eliptic wing formula
    C_l = 2*pi / (1 + 2 / AR) * deg2rad(alpha - alpha_low);
    Cd_i = 1 / (pi * AR) * C_l^2;
    alpha_i_rad = deg2rad(alpha - alpha_low) / (1 + AR/2);
end

function An = compute_An(n, b, m0, C, alpha, alpha_L0, theta)
    % SOLVE_COEFFICIENTS Solves for the first n coefficients A_n
    % Inputs:
    %   n        - Number of coefficients to compute
    %   b        - Given constant parameter
    %   m0       - Function handle or vector of m0(theta)
    %   C        - Function handle or vector of C(theta)
    %   alpha    - Function handle or vector of alpha(theta)
    %   alpha_L0 - Function handle or vector of alpha_L0(theta)
    %   theta    - Vector of theta values for numerical evaluation

    
    % Compute left-hand side scaling factor
    LHS_factor = -4 * b ./ (m0 .* C);
    
    % Construct the system of equations
    S = zeros(length(theta), n);
    T = zeros(length(theta), n);

    for k = 1:n
        S(:, k) = sin(k * theta);
        T(:, k) = k * (sin(k * theta) ./ sin(theta));
    end
    
    % Right-hand side of the equation
    RHS = zeros(length(theta),1);
    RHS(:) = -deg2rad(alpha - alpha_L0);
    
    % Solve the system
    A_matrix = (LHS_factor*S) - T;
    An = A_matrix \  RHS;

end