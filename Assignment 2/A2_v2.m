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


% fig1 = figure;
% set(fig1, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
% hold on
% yline(0, 'k', 'LineWidth', 2);
% for i = 1:length(Aspect_Ratio)
%     plot (AoA, C_l(i, :), LineWidth=3)
% end
% hold off
% grid on
% legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex')
% xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
% ylabel('Lift coefficient [-]', 'FontSize', 14, 'Interpreter', 'latex')
% title('Lift coefficient vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
% exportgraphics(fig1, 'Q1_Induced_CL_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');
% 
% fig2 = figure;
% set(fig2, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
% figure (2)
% hold on
% yline(0, 'k', 'LineWidth', 2);
% for i = 1:length(Aspect_Ratio)
%     plot (AoA, Cd_i(i, :), LineWidth=3)
% end
% hold off
% grid on
% legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
% xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
% ylabel('Induced Drag coefficent [-]', 'FontSize', 14, 'Interpreter', 'latex')
% title('Induced Drag coefficent vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
% exportgraphics(fig2, 'Q1_Induced_Cd_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');
% 
% fig3 = figure;
% set(fig3, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
% hold on
% yline(0, 'k', 'LineWidth', 2);
% for i = 1:length(Aspect_Ratio)
%     plot (AoA, CD(i, :), LineWidth=3)
% end
% hold off
% grid on
% legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
% xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
% ylabel('Drag coefficent [-]', 'FontSize', 14, 'Interpreter', 'latex')
% title('Drag coefficent vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
% exportgraphics(fig3, 'Q1_Cd_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');
% 
% fig4 = figure;
% set(fig4, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
% hold on
% yline(0, 'k', 'LineWidth', 2);
% for i = 1:length(Aspect_Ratio)
%     plot (AoA, alpha_i(i, :), LineWidth=3)
% end
% hold off
% grid on
% legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
% xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
% ylabel('Induced angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
% title('Induced angle of attack vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
% exportgraphics(fig4, 'Q1_Induced_AoA_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');



%% Question 2 Calculations

AoA_Q2 = [0, 5, 10];

m0 = 6.5958;
C = 1;

n = 50;
B = Aspect_Ratio .* C;
theta = linspace(0.01, pi, 180);
A_interm = zeros(n, length(AoA_Q2), length(B));


for i = 1:length(AoA_Q2)
    for j = 1:length(B)
        A_interm(:, i, j) = compute_An(n, B(j), m0, C, AoA_Q2(i), alpha_low, theta);
    end
end
Num = linspace(1,n,n);
An_plot1 = zeros(length(Num),length(theta));
for i = 1:length(theta)
    for j = 1:length(Num)
        An_plot1(j,i) = j*A_interm(j,2,1)*sin(Num(j)*theta(i));
    end
    
end
An_interm = sum(An_plot1, 1);

An_plot = 2*B(2)*1*An_interm;


alpha_i_q2 = zeros(length(AoA_Q2), length(B), length(theta));
for k = 1:length(theta)
    for i = 1:length(AoA_Q2)
        for j = 1:length(B)
            alpha_interm_sum = 0;  % Reset sum accumulator for each (i,j,k) combination
            for z = 1:n
                alpha_interm = Num(z) * A_interm(z,i,j) * (sin(Num(z)*theta(k)) / sin(theta(k)));
                alpha_interm_sum = alpha_interm_sum + alpha_interm;  % Add to the sum
            end
            alpha_i_q2(i,j,k) = alpha_interm_sum;  % Assign the final sum after the loop over z
            disp(['Theta:', num2str(rad2deg(theta(k))), ' | Alpha_i:', num2str(rad2deg(alpha_interm_sum))]);

        end
    end
end


fig5 = figure;
set(fig5, 'Units', 'pixels', 'Position', [100, 100, 1200, 900]);

for i = 1:length(AoA_Q2)
    subplot(length(AoA_Q2), 1, i);
    hold on;
    for j = 1:length(B)
        plot(rad2deg(theta), rad2deg(squeeze(alpha_i_q2(i, j, :))), 'DisplayName', ['AR = ' num2str(Aspect_Ratio(j))], LineWidth=2); 
    end
    title(['Indiced AoA vs. Theta for AoA = ' num2str(AoA_Q2(i)) ' deg'], 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('Theta (deg)', 'FontSize', 14, 'Interpreter', 'latex');
    grid on;
    ylabel('Induced AoA (deg)', 'FontSize', 14, 'Interpreter', 'latex');
    legend ('Location', 'north', 'FontSize', 12, 'Interpreter', 'latex');
    hold off;
end
sgtitle('Induced AoA vs. Theta for difefrent Angles of Attack', 'FontSize', 16, 'Interpreter', 'latex');
exportgraphics(fig5, 'Q2_Induced_AoA_vs_Theta.pdf', 'Resolution', 450, 'ContentType', 'vector');




%% For lift 

AoA_lift = linspace(-4, 10, 15);
A_1_lift = zeros(1, length(AoA_lift), length(B));
CL_Q2 = zeros(length(AoA_lift), length(B));


for i = 1:length(AoA_lift)
    for j = 1:length(B)
        A_1_lift(:, i, j) = compute_An(1, B(j), m0, C, AoA_lift(i), alpha_low, theta);
        CL_Q2(i,j) = pi*Aspect_Ratio(j)*A_1_lift(:,i,j);
    end
end

% fig6 = figure;
% set(fig6, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
% 
% hold on;
% for j = 1:length(B)
%     plot(AoA_lift, CL_Q2(:, j), 'DisplayName', ['AR = ' num2str(Aspect_Ratio(j))], LineWidth=3);
% end
% 
% xlabel('Angle of Attack (deg)', 'FontSize', 14, 'Interpreter', 'latex');
% ylabel('Coefficient of Lift (-)', 'FontSize', 14, 'Interpreter', 'latex');
% title('CL vs. AoA for Different Aspect Ratios', 'FontSize', 16, 'Interpreter', 'latex');
% legend ('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
% grid on;
% hold off;
% exportgraphics(fig6, 'Q2_CL_vs_AoA.pdf', 'Resolution', 500, 'ContentType', 'vector');


%% For Drag

n_drag = 50;
A_n_drag = zeros(n_drag, length(AoA_lift), length(B));
dell_Q2 = zeros(length(AoA_lift), length(B));
CD_Q2 = zeros(length(AoA_lift), length(B));

for i = 1:length(AoA_lift)
    for j = 1:length(B)
        A_n_drag(:, i, j) = compute_An(n_drag, B(j), m0, C, AoA_lift(i), alpha_low, theta);
        
        dell_sum = 0;  % Reset dell_sum for each iteration of i or j
        for k = 2:n_drag
            dell_const = k * ((A_n_drag(k,i,j) / A_n_drag(1,i,j))^2);
            dell_sum = dell_sum + dell_const;
        end
        
        dell_Q2(i,j) = dell_sum;
        CD_Q2(i,j) = ((CL_Q2(i,j)^2) / (pi * Aspect_Ratio(j))) * (1 + dell_Q2(i,j));
    end
end

% fig7 = figure;
% set(fig7, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
% hold on;
% for j = 1:length(B)
%     plot(AoA_lift, CD_Q2(:, j), 'DisplayName', ['AR = ' num2str(Aspect_Ratio(j))], LineWidth=3);  % Plot CD_Q2 for each B
% end
% 
% xlabel('Angle of Attack (deg)', 'FontSize', 14, 'Interpreter', 'latex');
% ylabel('Coefficient of Drag (-)', 'FontSize', 14, 'Interpreter', 'latex');
% title('CD vs. AoA for Different Aspect Ratios', 'FontSize', 16, 'Interpreter', 'latex');
% grid on;
% legend ('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
% hold off;
% exportgraphics(fig7, 'Q2_CD_vs_AoA.pdf', 'Resolution', 500, 'ContentType', 'vector');


close all   % DELETE THIS COMMAND =============================================================================================


%% Question 3


% Given parameters
AR = 6;                 % Aspect Ratio
b = 1;                  % Assume normalized wingspan b = 1 for simplicity
TR_values = [0.2, 0.4, 0.6, 0.8, 1.0]; % Taper ratios
AOA = 60;               % Angle of attack in degrees
U_inf = 1;              % Normalized free-stream velocity
N = 100;                % Number of spanwise points
x_tilda = linspace(-1, 1, N); % Non-dimensional spanwise coordinate

% Loop over different taper ratios
fig8 = figure;
set(fig8, 'Units', 'pixels', 'Position', [100, 100, 1200, 900]);

for TR = TR_values
    % Compute chord distribution
    c_tilda = (2 * (1 - (1 - TR) * abs(x_tilda))) / (1 + TR);

    % Compute circulation distribution (simplified as a function of chord)
    Gamma_tilda = c_tilda .* sind(AOA); % Proportional to local lift
    
    % Compute induced angle of attack (simplified relation)
    alpha_i = Gamma_tilda / (pi * b/2); % Approximate induced AOA in radians
    alpha_i_deg = rad2deg(alpha_i);     % Convert to degrees
    
    % Compute lift coefficient
    C_l = 2 * Gamma_tilda ./ c_tilda;
    
    % Compute induced drag coefficient
    CD_i = C_l .* sind(alpha_i_deg); % Approximate relation for induced drag

    % Plot results
    subplot(2,2,1);
    plot(x_tilda, Gamma_tilda, 'DisplayName', sprintf('TR = %.1f', TR), LineWidth=2);
    hold on;
    xlabel('$\tilde{x}$', 'Interpreter', 'latex'); ylabel('$\tilde{\Gamma}$', 'Interpreter', 'latex');
    title('Dimensionless Circulation Distribution');
    grid on;
    legend ('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');

    subplot(2,2,2);
    plot(x_tilda, alpha_i_deg, 'DisplayName', sprintf('TR = %.1f', TR), LineWidth=2);
    hold on;
    xlabel('$\tilde{x}$', 'Interpreter', 'latex'); ylabel('$\alpha_i$ [deg]', 'Interpreter', 'latex');
    title('Induced Angle of Attack');
    grid on;
    legend ('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');

    subplot(2,2,3);
    plot(x_tilda, C_l, 'DisplayName', sprintf('TR = %.1f', TR), LineWidth=2);
    hold on;
    xlabel('$\tilde{x}$', 'Interpreter', 'latex'); ylabel('$C_l$', 'Interpreter', 'latex');
    title('Lift Coefficient Distribution');
    grid on;
    legend ('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');

    subplot(2,2,4);
    plot(x_tilda, CD_i, 'DisplayName', sprintf('TR = %.1f', TR), LineWidth=2);
    hold on;
    xlabel('$\tilde{x}$', 'Interpreter', 'latex'); ylabel('$C_{D_i}$', 'Interpreter', 'latex');
    title('Induced Drag Coefficient Distribution');
    grid on;
    legend  ('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
end
sgtitle('Parameters along the spanwise coordinates', 'FontSize', 16, 'Interpreter', 'latex');
exportgraphics(fig8, 'Q3_Induced_AoA_vs_Theta.pdf', 'Resolution', 450, 'ContentType', 'vector');



%% Functions

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