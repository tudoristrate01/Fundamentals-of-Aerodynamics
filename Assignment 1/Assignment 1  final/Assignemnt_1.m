clc 
clear all
close all

%% Question 1

t_2312 = 0.12;  % relative thickness
m_2312 = 0.02;  % max camber
p_2312 = 0.3;  % location of the max camber from leading edge

t_2324 = 0.24;
m_2324 = 0.02;
p_2324 = 0.3;

t_4412 = 0.12;
m_4412 = 0.04;
p_4412 = 0.4;

t_4424 = 0.24;
m_4424 = 0.04;
p_4424 = 0.4;

% Compute airfoil geometries
[x_c_2312, y_c_2312, x_u_2312, y_u_2312, x_l_2312, y_l_2312, theta_2312] = plot_airfoil(t_2312, m_2312, p_2312);
[x_c_2324, y_c_2324, x_u_2324, y_u_2324, x_l_2324, y_l_2324, theta_2324] = plot_airfoil(t_2324, m_2324, p_2324);
[x_c_4412, y_c_4412, x_u_4412, y_u_4412, x_l_4412, y_l_4412, theta_4412] = plot_airfoil(t_4412, m_4412, p_4412);
[x_c_4424, y_c_4424, x_u_4424, y_u_4424, x_l_4424, y_l_4424, theta_4424] = plot_airfoil(t_4424, m_4424, p_4424);

figure(1)
hold on
plot(x_l_2312, y_l_2312, 'Color', 'b');
h1 = plot(x_u_2312, y_u_2312, 'Color', 'b', 'DisplayName', 'NACA 2312');
plot(x_c_2312, y_c_2312, 'Color', 'b', LineStyle='--',DisplayName='Camber 2312' )
plot(x_l_2324, y_l_2324, 'Color', 'r');
h2 = plot(x_u_2324, y_u_2324, 'Color', 'r', 'DisplayName', 'NACA 2324');
plot(x_c_2324, y_c_2324, 'Color', 'r', LineStyle='-.', DisplayName='Camber 2324')
hold off
grid on
axis equal
legend([h1, h2])
xlabel('Chordwise Position')
ylabel('Thickness / Camber')
title('NACA 4412 vs NACA 4424')

figure(2)
hold on
plot (x_l_4412, y_l_4412, 'b');
h3 = plot(x_u_4412, y_u_4412, 'b', DisplayName='NACA 4412');
plot(x_c_4412, y_c_4412, '--b')
plot (x_l_4424, y_l_4424, 'r');
h4 = plot(x_u_4424, y_u_4424, 'r', DisplayName='NACA 4424');
plot(x_c_4424, y_c_4424, '-.r')
hold off
grid on
axis equal
legend([h3, h4])
xlabel('Chordwise Position')
ylabel('Thickness / Camber')
title('NACA 4412 vs NACA 4424')

%% Question 2 - Thin Airfoil Model

% Thin airfoil theory
alpha = linspace(-10,15,100);  % Creates 100 points between -10 and 15

% NACA 2312
[dy_dxc_theta_2312, A0_2312, A1_2312, c_lift_2312] = thin_airfoil_deta_dx(p_2312, alpha, m_2312);


% NACA 2324
[dy_dxc_theta_2324, A0_2324, A1_2324, c_lift_2324] = thin_airfoil_deta_dx(p_2324, alpha, m_2324);



% NACA 4412
[dy_dxc_theta_4412, A0_4412, A1_4412, c_lift_4412] = thin_airfoil_deta_dx(p_4412, alpha, m_4412);

% NACA 4424
[dy_dxc_theta_4424, A0_4424, A1_4424, c_lift_4424] = thin_airfoil_deta_dx(p_4424, alpha, m_4424);


%% Question 2 - Panel Method

% NACA 2312

[x_c_panel_2312, y_c_panel_2312] = orient_airfoil(x_u_2312, y_u_2312, x_l_2312, y_l_2312);
Cl_2312_panel = zeros(length(alpha),1);

for i = 1:length(alpha)
    [Gamma_new_2312, Cl_2312_old, CP_2312_q2, x_cp_2312] = Panel_method_new(alpha(i), 0, x_c_panel_2312, y_c_panel_2312,1);
    [Gamma_2312, Cl_2312_panel(i), CP_2312_q2, x_cp_2312] = Panel_method_new(alpha(i), Gamma_new_2312, x_c_panel_2312, y_c_panel_2312,1);
end

% NACA 2324

[x_c_panel_2324, y_c_panel_2324] = orient_airfoil(x_u_2324, y_u_2324, x_l_2324, y_l_2324);
Cl_2324_panel = zeros(length(alpha),1);

for i = 1:length(alpha)
    [Gamma_new_2324, Cl_2324_old, CP_2324_q2, x_cp_2324] = Panel_method_new(alpha(i), 0, x_c_panel_2324, y_c_panel_2324,1);
    [Gamma_2324, Cl_2324_panel(i), CP_2324_q2, x_cp_2324] = Panel_method_new(alpha(i), Gamma_new_2324, x_c_panel_2324, y_c_panel_2324,1);
end

% NACA 4412

[x_c_panel_4412, y_c_panel_4412] = orient_airfoil(x_u_4412, y_u_4412, x_l_4412, y_l_4412);
Cl_4412_panel = zeros(length(alpha),1);

for i = 1:length(alpha)
    [Gamma_new_4412, Cl_4412_old, CP_4412_q2, x_cp_4412] = Panel_method_new(alpha(i), 0, x_c_panel_4412, y_c_panel_4412,1);
    [Gamma_4412, Cl_4412_panel(i), CP_4412_q2, x_cp_4412] = Panel_method_new(alpha(i), Gamma_new_4412, x_c_panel_4412, y_c_panel_4412,1);
end

% NACA 4424

[x_c_panel_4424, y_c_panel_4424] = orient_airfoil(x_u_4424, y_u_4424, x_l_4424, y_l_4424);
Cl_4424_panel = zeros(length(alpha),1);

for i = 1:length(alpha)
    [Gamma_new_4424, Cl_4424_old, CP_4424_q2, x_cp_4424] = Panel_method_new(alpha(i), 0, x_c_panel_4424, y_c_panel_4424,1);
    [Gamma_4424, Cl_4424_panel(i), CP_4424_q2, x_cp_4424] = Panel_method_new(alpha(i), Gamma_new_4424, x_c_panel_4424, y_c_panel_4412,1);
end


%% Question 2 - XFOIL (Free)

% NACA 2312
% Read the data from the text file
data_2312_free = load('NACA_2312_free.txt');  

% Assign the first and second columns to the variables
alpha_2312_xfoil_free = data_2312_free(:, 1);   % First column for x
cl_2312_xfoil_free = data_2312_free(:, 2);  % Second column for cl

% NACA 2324
% Read the data from the text file
data_2324_free = load('NACA_2324_free.txt');  

% Assign the first and second columns to the variables
alpha_2324_xfoil_free = data_2324_free(:, 1);   % First column for x
cl_2324_xfoil_free = data_2324_free(:, 2);  % Second column for cl

% NACA 4412
% Read the data from the text file
data_4412_free = load('NACA_4412_free.txt');  

% Assign the first and second columns to the variables
alpha_4412_xfoil_free = data_4412_free(:, 1);   % First column for x
cl_4412_xfoil_free = data_4412_free(:, 2);  % Second column for cl

% NACA 2312
% Read the data from the text file
data_4424_free = load('NACA_4424_free.txt');  

% Assign the first and second columns to the variables
alpha_4424_xfoil_free = data_4424_free(:, 1);   % First column for x
cl_4424_xfoil_free = data_4424_free(:, 2);  % Second column for cl

%% Question 2 - XFOIL (Fixed)

% NACA 2312
% Read the data from the text file
data_2312_fixed = load('NACA_2312_fixed.txt');  

% Assign the first and second columns to the variables
alpha_2312_xfoil_fixed = data_2312_fixed(:, 1);   % First column for x
cl_2312_xfoil_fixed = data_2312_fixed(:, 2);  % Second column for cl

% NACA 2324
% Read the data from the text file
data_2324_fixed = load('NACA_2324_fixed.txt');  

% Assign the first and second columns to the variables
alpha_2324_xfoil_fixed = data_2324_fixed(:, 1);   % First column for x
cl_2324_xfoil_fixed = data_2324_fixed(:, 2);  % Second column for cl

% NACA 4412
% Read the data from the text file
data_4412_fixed = load('NACA_4412_fixed.txt');  

% Assign the first and second columns to the variables
alpha_4412_xfoil_fixed = data_4412_fixed(:, 1);   % First column for x
cl_4412_xfoil_fixed = data_4412_fixed(:, 2);  % Second column for cl

% NACA 2312
% Read the data from the text file
data_4424_fixed = load('NACA_4424_fixed.txt');  

% Assign the first and second columns to the variables
alpha_4424_xfoil_fixed = data_4424_fixed(:, 1);   % First column for x
cl_4424_xfoil_fixed = data_4424_fixed(:, 2);  % Second column for cl


%% Question 2 - Plots

figure;

% Define airfoil names for titles
airfoil_names = {'NACA 2312', 'NACA 2324', 'NACA 4412', 'NACA 4424'};

data_sets = {...
    {alpha, c_lift_2312, Cl_2312_panel, alpha_2312_xfoil_free, cl_2312_xfoil_free, alpha_2312_xfoil_fixed, cl_2312_xfoil_fixed},...
    {alpha, c_lift_2324, Cl_2324_panel, alpha_2324_xfoil_free, cl_2324_xfoil_free, alpha_2324_xfoil_fixed, cl_2324_xfoil_fixed},...
    {alpha, c_lift_4412, Cl_4412_panel, alpha_4412_xfoil_free, cl_4412_xfoil_free, alpha_4412_xfoil_fixed, cl_4412_xfoil_fixed},...
    {alpha, c_lift_4424, Cl_4424_panel, alpha_4424_xfoil_free, cl_4424_xfoil_free, alpha_4424_xfoil_fixed, cl_4424_xfoil_fixed}...
};

% Loop over each airfoil and create subplots
for i = 1:4
    subplot(2,2,i);
    hold on; grid on; box on;
    
    % Extract data for the current airfoil
    alpha_data = data_sets{i}{1};
    c_lift = data_sets{i}{2};
    Cl_panel = data_sets{i}{3};
    alpha_xfoil_free = data_sets{i}{4};
    cl_xfoil_free = data_sets{i}{5};
    alpha_xfoil_fixed = data_sets{i}{6};
    cl_xfoil_fixed = data_sets{i}{7};
    
    % Plot data
    plot(alpha_data, c_lift, 'b-', 'LineWidth', 2);
    plot(alpha_data, Cl_panel, 'g-', 'LineWidth', 2);
    plot(alpha_xfoil_free, cl_xfoil_free, 'r-', 'LineWidth', 2);
    plot(alpha_xfoil_fixed, cl_xfoil_fixed, 'm-', 'LineWidth', 2);
    
    % Labels and title
    xlabel('Angle of Attack, \alpha (deg)', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('Lift Coefficient, C_L', 'FontSize', 12, 'Interpreter', 'latex');
    title(['Lift Coefficient vs. Angle of Attack for ', airfoil_names{i}], 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    % Add X and Y Axes
    xline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Customize grid and axis
    ax = gca;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.FontSize = 12;
    ax.LineWidth = 1.2;
end

% Add a single legend for all subplots
legend({'Thin Airfoil Theory', 'Panel Method', 'XFOIL - Free', 'XFOIL - Fixed'}, 'Location', 'BestOutside', 'FontSize', 12, 'Interpreter', 'latex', 'EdgeColor', 'k');

%% Question 3 - Thin Airfoil theory

rho = 1.225;
U0 = 1;
alpha_q3 = 10;
N = 50;

% NACA 2312
[x_cp_2312, dellcP_q3_2312] = thin_airfoil_dellcP(dy_dxc_theta_2312, alpha_q3, U0, N);

% NACA 2324
[x_cp_2324, dellcP_q3_2324] = thin_airfoil_dellcP(dy_dxc_theta_2324, alpha_q3, U0, N);

% NACA 4412
[x_cp_4412, dellcP_q3_4412] = thin_airfoil_dellcP(dy_dxc_theta_4412, alpha_q3, U0, N);

% NACA 4424
[x_cp_4424, dellcP_q3_4424] = thin_airfoil_dellcP(dy_dxc_theta_4424, alpha_q3, U0, N);

%% Question 3 - Panel method

% NACA 2312
[x_u_2312_new, y_u_2312_new, x_l_2312_new, y_l_2312_new] = sort_x_y_panel(x_u_2312, y_u_2312, x_l_2312, y_l_2312);

[Gamma_new_2312_u_q3, Cl_2312_old_u_q3, CP_2312_old_u_q3, x_cp_2312_old_u_q3] = Panel_method_new(10, 0, x_u_2312_new', y_u_2312_new',1);
[Gamma_2312_u_q3, Cl_2312_panel_u_q3, CP_2312_u_q3, x_cp_2312_u_q3] = Panel_method_new(10, Gamma_new_2312_u_q3, x_u_2312_new', y_u_2312_new',1);

[Gamma_new_2312_l_q3, Cl_2312_old_l_q3, CP_2312_old_l_q3, x_cp_2312_old_l_q3] = Panel_method_new(10, 0, x_l_2312_new', y_l_2312_new',1);
[Gamma_2312_l_q3, Cl_2312_panel_l_q3, CP_2312_l_q3, x_cp_2312_l_q3] = Panel_method_new(10, Gamma_new_2312_l_q3, x_l_2312_new', y_l_2312_new',1);

[x_cp_panel_2312, dell_CP_panel_2312] = cP_Panel(x_cp_2312_l_q3, CP_2312_l_q3, x_cp_2312_u_q3, CP_2312_u_q3);

% NACA 2324
[x_u_2324_new, y_u_2324_new, x_l_2324_new, y_l_2324_new] = sort_x_y_panel(x_u_2324, y_u_2324, x_l_2324, y_l_2324);

[Gamma_new_2324_u_q3, Cl_2324_old_u_q3, CP_2324_old_u_q3, x_cp_2324_old_u_q3] = Panel_method_new(10, 0, x_u_2324_new', y_u_2324_new',1);
[Gamma_2324_u_q3, Cl_2324_panel_u_q3, CP_2324_u_q3, x_cp_2324_u_q3] = Panel_method_new(10, Gamma_new_2324_u_q3, x_u_2324_new', y_u_2324_new',1);

[Gamma_new_2324_l_q3, Cl_2324_old_l_q3, CP_2324_old_l_q3, x_cp_2324_old_l_q3] = Panel_method_new(10, 0, x_l_2324_new', y_l_2324_new',1);
[Gamma_2324_l_q3, Cl_2324_panel_l_q3, CP_2324_l_q3, x_cp_2324_l_q3] = Panel_method_new(10, Gamma_new_2324_l_q3, x_l_2324_new', y_l_2324_new',1);

[x_cp_panel_2324, dell_CP_panel_2324] = cP_Panel(x_cp_2324_l_q3, CP_2324_l_q3, x_cp_2324_u_q3, CP_2324_u_q3);

% NACA 4412
[x_u_4412_new, y_u_4412_new, x_l_4412_new, y_l_4412_new] = sort_x_y_panel(x_u_4412, y_u_4412, x_l_4412, y_l_4412);

[Gamma_new_4412_u_q3, Cl_4412_old_u_q3, CP_4412_old_u_q3, x_cp_4412_old_u_q3] = Panel_method_new(10, 0, x_u_4412_new', y_u_4412_new',1);
[Gamma_4412_u_q3, Cl_4412_panel_u_q3, CP_4412_u_q3, x_cp_4412_u_q3] = Panel_method_new(10, Gamma_new_4412_u_q3, x_u_4412_new', y_u_4412_new',1);

[Gamma_new_4412_l_q3, Cl_4412_old_l_q3, CP_4412_old_l_q3, x_cp_4412_old_l_q3] = Panel_method_new(10, 0, x_l_4412_new', y_l_4412_new',1);
[Gamma_4412_l_q3, Cl_4412_panel_l_q3, CP_4412_l_q3, x_cp_4412_l_q3] = Panel_method_new(10, Gamma_new_4412_l_q3, x_l_4412_new', y_l_4412_new',1);

[x_cp_panel_4412, dell_CP_panel_4412] = cP_Panel(x_cp_4412_l_q3, CP_4412_l_q3, x_cp_4412_u_q3, CP_4412_u_q3);


% NACA 4424
[x_u_4424_new, y_u_4424_new, x_l_4424_new, y_l_4424_new] = sort_x_y_panel(x_u_4424, y_u_4424, x_l_4424, y_l_4424);

[Gamma_new_4424_u_q3, Cl_4424_old_u_q3, CP_4424_old_u_q3, x_cp_4424_old_u_q3] = Panel_method_new(10, 0, x_u_4424_new', y_u_4424_new',1);
[Gamma_4424_u_q3, Cl_4424_panel_u_q3, CP_4424_u_q3, x_cp_4424_u_q3] = Panel_method_new(10, Gamma_new_4424_u_q3, x_u_4424_new', y_u_4424_new',1);

[Gamma_new_4424_l_q3, Cl_4424_old_l_q3, CP_4424_old_l_q3, x_cp_4424_old_l_q3] = Panel_method_new(10, 0, x_l_4424_new', y_l_4424_new',1);
[Gamma_4424_l_q3, Cl_4424_panel_l_q3, CP_4424_l_q3, x_cp_4424_l_q3] = Panel_method_new(10, Gamma_new_4424_l_q3, x_l_4424_new', y_l_4424_new',1);

[x_cp_panel_4424, dell_CP_panel_4424] = cP_Panel(x_cp_4424_l_q3, CP_4424_l_q3, x_cp_4424_u_q3, CP_4424_u_q3);


%% Question 3 - XFOIL

% Read the data from the text file
data_2312_free_lower_CP = load('NACA_2312_lower_CP.txt');  

% Assign the first and second columns to the variables
x_lower_2312_CP = data_2312_free_lower_CP(:, 1);   % First column for x
y_lower_2312_CP = data_2312_free_lower_CP(:, 2);  % Second column for cl
CP_lower_2312_CP = data_2312_free_lower_CP(:, 3);


% Read the data from the text file
data_2312_free_upper_CP = load('NACA_2312_upper_CP.txt');  

% Assign the first and second columns to the variables
x_upper_2312_CP = data_2312_free_upper_CP(:, 1);   % First column for x
y_upper_2312_CP = data_2312_free_upper_CP(:, 2);  % Second column for cl
CP_upper_2312_CP = data_2312_free_upper_CP(:, 3);

[x_u_sorted_2312, y_u_sorted_2312, CP_u_sorted_2312, x_l_sorted_2312, y_l_sorted_2312, CP_l_sorted_2312] = sort_xfoil_CP(x_upper_2312_CP, y_upper_2312_CP, CP_upper_2312_CP, x_lower_2312_CP, y_lower_2312_CP, CP_lower_2312_CP);
[x_cp_xfoil_2312, dell_CP_xfoil_2312] = cP_Panel(x_l_sorted_2312, CP_l_sorted_2312, x_u_sorted_2312, CP_u_sorted_2312);


%% Question 3 - Plots

% NACA 2312

figure;
% Plot Thin Airfoil Theory with markers and thicker line
plot(x_cp_2312(2:end-1), dellcP_q3_2312(2:end-1), 'b-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% Plot Panel Method with markers and thicker line
plot(x_cp_panel_2312(2:end-1), dell_CP_panel_2312(2:end-1), 'r-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(x_cp_xfoil_2312(2:end-1), dell_CP_xfoil_2312(2:end-1), 'g-', 'LineWidth', 1.5, 'MarkerSize', 6);


% Adding grid for better readability
grid on;
% Label the x and y axes
xlabel('x/c', 'FontSize', 14);
ylabel('$\Delta C_P$', 'Interpreter', 'latex', 'FontSize', 14);
% Add a title with latex formatting
title('$\Delta C_P$ vs. $x/c$ - NACA 2312', 'Interpreter', 'latex', 'FontSize', 16);
% Add a legend with a background and adjusted font size
legend('Thin airfoil theory', 'Panel method', 'XFOIL-free','Location', 'best', 'FontSize', 12, 'Box', 'on');
hold off;


% NACA 2324

figure;
% Plot Thin Airfoil Theory with markers and thicker line
plot(x_cp_2324(2:end-1), dellcP_q3_2324(2:end-1), 'b-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% Plot Panel Method with markers and thicker line
plot(x_cp_panel_2324(2:end-1), dell_CP_panel_2324(2:end-1), 'r-', 'LineWidth', 1.5, 'MarkerSize', 6);
% Adding grid for better readability
grid on;
% Label the x and y axes
xlabel('x/c', 'FontSize', 14);
ylabel('$\Delta C_P$', 'Interpreter', 'latex', 'FontSize', 14);
% Add a title with latex formatting
title('$\Delta C_P$ vs. $x/c$ - NACA 2324', 'Interpreter', 'latex', 'FontSize', 16);
% Add a legend with a background and adjusted font size
legend('Thin airfoil theory', 'Panel method', 'Location', 'best', 'FontSize', 12, 'Box', 'on');
hold off;


% NACA 4412

figure;
% Plot Thin Airfoil Theory with markers and thicker line
plot(x_cp_4412(2:end-1), dellcP_q3_4412(2:end-1), 'b-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% Plot Panel Method with markers and thicker line
plot(x_cp_panel_4412(2:end-1), dell_CP_panel_4412(2:end-1), 'r-', 'LineWidth', 1.5, 'MarkerSize', 6);
% Adding grid for better readability
grid on;
% Label the x and y axes
xlabel('x/c', 'FontSize', 14);
ylabel('$\Delta C_P$', 'Interpreter', 'latex', 'FontSize', 14);
% Add a title with latex formatting
title('$\Delta C_P$ vs. $x/c$ - NACA 4412', 'Interpreter', 'latex', 'FontSize', 16);
% Add a legend with a background and adjusted font size
legend('Thin airfoil theory', 'Panel method', 'Location', 'best', 'FontSize', 12, 'Box', 'on');
hold off;


% NACA 4424

figure;
% Plot Thin Airfoil Theory with markers and thicker line
plot(x_cp_4424(2:end-1), dellcP_q3_4424(2:end-1), 'b-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% Plot Panel Method with markers and thicker line
plot(x_cp_panel_4424(2:end-1), dell_CP_panel_4424(2:end-1), 'r-', 'LineWidth', 1.5, 'MarkerSize', 6);
% Adding grid for better readability
grid on;
% Label the x and y axes
xlabel('x/c', 'FontSize', 14);
ylabel('$\Delta C_P$', 'Interpreter', 'latex', 'FontSize', 14);
% Add a title with latex formatting
title('$\Delta C_P$ vs. $x/c$ - NACA 4424', 'Interpreter', 'latex', 'FontSize', 16);
% Add a legend with a background and adjusted font size
legend('Thin airfoil theory', 'Panel method', 'Location', 'best', 'FontSize', 12, 'Box', 'on');
hold off;

%% Question 4 - Panel method

% NACA 2312
[Gamma_new_2312_q4, Cl_2312_old_q4, CP_2312_old_q4, x_cp_2312_old_q4] = Panel_method_new(10, 0, x_c_panel_2312, y_c_panel_2312,1);
[Gamma_2312_q4, Cl_2312_panel_q4, CP_2312_q4, x_cp_2312_q4] = Panel_method_new(10, Gamma_new_2312_q4, x_c_panel_2312, y_c_panel_2312,1);

% NACA 2324
[Gamma_new_2324_q4, Cl_2324_old_q4, CP_2324_old_q4, x_cp_2324_old_q4] = Panel_method_new(10, 0, x_c_panel_2324, y_c_panel_2324,1);
[Gamma_2324_q4, Cl_2324_panel_q4, CP_2324_q4, x_cp_2324_q4] = Panel_method_new(10, Gamma_new_2324_q4, x_c_panel_2324, y_c_panel_2324,1);

% NACA 4412
[Gamma_new_4412_q4, Cl_4412_old_q4, CP_4412_old_q4, x_cp_4412_old_q4] = Panel_method_new(10, 0, x_c_panel_4412, y_c_panel_4412,1);
[Gamma_4412_q4, Cl_4412_panel_q4, CP_4412_q4, x_cp_4412_q4] = Panel_method_new(10, Gamma_new_4412_q4, x_c_panel_4412, y_c_panel_4412,1);

% NACA 4424
[Gamma_new_4424_q4, Cl_4424_old_q4, CP_4424_old_q4, x_cp_4424_old_q4] = Panel_method_new(10, 0, x_c_panel_4424, y_c_panel_4424,1);
[Gamma_4424_q4, Cl_4424_panel_q4, CP_4424_q4, x_cp_4424_q4] = Panel_method_new(10, Gamma_new_4424_q4, x_c_panel_4424, y_c_panel_4424,1);

%% Question 4 - XFOIL (Free)

% NACA 2312
% Read the data from the text file
data_2312_free_CP = load('NACA_2312_free_CP.txt');  

% Assign the first and second columns to the variables
x_2312_xfoil_free = data_2312_free_CP(:, 1);   % First column for x
cP_2312_xfoil_free = data_2312_free_CP(:, 3);  

% NACA 2324
% Read the data from the text file
data_2324_free_CP = load('NACA_2324_free_CP.txt');  

% Assign the first and second columns to the variables
x_2324_xfoil_free = data_2324_free_CP(:, 1);   % First column for x
cP_2324_xfoil_free = data_2324_free_CP(:, 3);  

% NACA 4412
% Read the data from the text file
data_4412_free_CP = load('NACA_4412_free_CP.txt');  

% Assign the first and second columns to the variables
x_4412_xfoil_free = data_4412_free_CP(:, 1);   % First column for x
cP_4412_xfoil_free = data_4412_free_CP(:, 3);  

% NACA 2312
% Read the data from the text file
data_4424_free_CP = load('NACA_4424_free_CP.txt');  

% Assign the first and second columns to the variables
x_4424_xfoil_free = data_4424_free_CP(:, 1);   % First column for x
cP_4424_xfoil_free = data_4424_free_CP(:, 3);  

%% Question 4 - XFOIL (Fixed)

% NACA 2312
% Read the data from the text file
data_2312_fixed_CP = load('NACA_2312_fixed_CP.txt');  

% Assign the first and second columns to the variables
x_2312_xfoil_fixed = data_2312_fixed_CP(:, 1);   % First column for x
cP_2312_xfoil_fixed = data_2312_fixed_CP(:, 3);  

% NACA 2324
% Read the data from the text file
data_2324_fixed_CP = load('NACA_2324_fixed_CP.txt');  

% Assign the first and second columns to the variables
x_2324_xfoil_fixed = data_2324_fixed_CP(:, 1);   % First column for x
cP_2324_xfoil_fixed = data_2324_fixed_CP(:, 3);  

% NACA 4412
% Read the data from the text file
data_4412_fixed_CP = load('NACA_4412_fixed_CP.txt');  

% Assign the first and second columns to the variables
x_4412_xfoil_fixed = data_4412_fixed_CP(:, 1);   % First column for x
cP_4412_xfoil_fixed = data_4412_fixed_CP(:, 3); 

% NACA 2312
% Read the data from the text file
data_4424_fixed_CP = load('NACA_4424_fixed_CP.txt');  

% Assign the first and second columns to the variables
x_4424_xfoil_fixed = data_4424_fixed_CP(:, 1);   % First column for x
cP_4424_xfoil_fixed = data_4424_fixed_CP(:, 3);  


%% Question 4 - Plots

% NACA 2312
% NACA 2312 - Pressure Coefficient Plot
figure;
hold on; grid on; box on;

% Plot Panel Method (Blue)
plot(x_cp_2312_q4, CP_2312_q4, 'b-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');

% Plot XFOIL - Free Transition (Red)
plot(x_2312_xfoil_free, cP_2312_xfoil_free, 'r-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');

% Plot XFOIL - Fixed Transition (Magenta)
plot(x_2312_xfoil_fixed, cP_2312_xfoil_fixed, 'k--', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm');

% Labels and Title with LaTeX formatting
xlabel('$x/c$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$C_P$ [-]', 'FontSize', 16, 'Interpreter', 'latex');
title('$C_P$ vs. $x/c$ for NACA 2312', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Invert Y-axis (Aerodynamic Standard)
set(gca, 'YDir', 'reverse');

% Add a Legend with Distinct Labels
legend({'Panel Method', 'XFOIL - Free Transition', 'XFOIL - Fixed Transition'}, ...
       'Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', ...
       'Box', 'on', 'EdgeColor', 'k');

% Customize grid and axis
ax = gca;
ax.XMinorGrid = 'on'; % Enable minor grid for better readability
ax.YMinorGrid = 'on';
ax.FontSize = 14; % Increase font size for better visibility
ax.LineWidth = 1.5; % Thicker axis lines

% Display plot
hold off;


% NACA 2324
% NACA 2324 - Pressure Coefficient Plot
figure;
hold on; grid on; box on;

% Plot Panel Method (Blue)
plot(x_cp_2324_q4, CP_2324_q4, 'b-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');

% Plot XFOIL - Free Transition (Red)
plot(x_2324_xfoil_free, cP_2324_xfoil_free, 'r-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');

% Plot XFOIL - Fixed Transition (Magenta)
plot(x_2324_xfoil_fixed, cP_2324_xfoil_fixed, 'k--', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm');

% Labels and Title with LaTeX formatting
xlabel('$x/c$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$C_P$ [-]', 'FontSize', 16, 'Interpreter', 'latex');
title('$C_P$ vs. $x/c$ for NACA 2324', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Invert Y-axis (Aerodynamic Standard)
set(gca, 'YDir', 'reverse');

% Add a Legend with Distinct Labels
legend({'Panel Method', 'XFOIL - Free Transition', 'XFOIL - Fixed Transition'}, ...
       'Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', ...
       'Box', 'on', 'EdgeColor', 'k');

% Customize grid and axis
ax = gca;
ax.XMinorGrid = 'on'; % Enable minor grid for better readability
ax.YMinorGrid = 'on';
ax.FontSize = 14; % Increase font size for better visibility
ax.LineWidth = 1.5; % Thicker axis lines

% Display plot
hold off;


% NACA 4412
figure;
hold on; grid on; box on;

% Plot Panel Method (Blue)
plot(x_cp_4412_q4, CP_4412_q4, 'b-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');

% Plot XFOIL - Free Transition (Red)
plot(x_4412_xfoil_free, cP_4412_xfoil_free, 'r-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');

% Plot XFOIL - Fixed Transition (Magenta)
plot(x_4412_xfoil_fixed, cP_4412_xfoil_fixed, 'k--', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm');

% Labels and Title with LaTeX formatting
xlabel('$x/c$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$C_P$ [-]', 'FontSize', 16, 'Interpreter', 'latex');
title('$C_P$ vs. $x/c$ for NACA 4412', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Invert Y-axis (Aerodynamic Standard)
set(gca, 'YDir', 'reverse');

% Add a Legend with Distinct Labels
legend({'Panel Method', 'XFOIL - Free Transition', 'XFOIL - Fixed Transition'}, ...
       'Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', ...
       'Box', 'on', 'EdgeColor', 'k');

% Customize grid and axis
ax = gca;
ax.XMinorGrid = 'on'; % Enable minor grid for better readability
ax.YMinorGrid = 'on';
ax.FontSize = 14; % Increase font size for better visibility
ax.LineWidth = 1.5; % Thicker axis lines

% Display plot
hold off;

% NACA 4424
figure;
hold on; grid on; box on;

% Plot Panel Method (Blue)
plot(x_cp_4424_q4, CP_4424_q4, 'b-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');

% Plot XFOIL - Free Transition (Red)
plot(x_4424_xfoil_free, cP_4424_xfoil_free, 'r-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');

% Plot XFOIL - Fixed Transition (Magenta)
plot(x_4424_xfoil_fixed, cP_4424_xfoil_fixed, 'k--', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm');

% Labels and Title with LaTeX formatting
xlabel('$x/c$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$C_P$ [-]', 'FontSize', 16, 'Interpreter', 'latex');
title('$C_P$ vs. $x/c$ for NACA 4424', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Invert Y-axis (Aerodynamic Standard)
set(gca, 'YDir', 'reverse');

% Add a Legend with Distinct Labels
legend({'Panel Method', 'XFOIL - Free Transition', 'XFOIL - Fixed Transition'}, ...
       'Location', 'NorthEast', 'FontSize', 14, 'Interpreter', 'latex', ...
       'Box', 'on', 'EdgeColor', 'k');

% Customize grid and axis
ax = gca;
ax.XMinorGrid = 'on'; % Enable minor grid for better readability
ax.YMinorGrid = 'on';
ax.FontSize = 14; % Increase font size for better visibility
ax.LineWidth = 1.5; % Thicker axis lines

% Display plot
hold off;

%% Functions

function [yc] = camber_line(m, p, x_c)
    if x_c <= p
        yc = m ./ p^2 .* (2 .* p .* x_c - x_c.^2);
    else
        yc = m ./ (1 - p).^2 .* (1 - 2.*p + 2 .* p .* x_c - x_c.^2);
    end
end

function [yt] = half_thickness(t, x)
    yt = 5*t * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x.^2 + 0.2843 * x.^3 - 0.1015 * x.^4);
end

function [xu, yu] = upper_line(x, yt, yc, theta)
    xu = x - yt .* sin(theta);
    yu = yc + yt .* cos(theta);
end

function [xl, yl] = lower_line(x, yt, yc, theta)
    xl = x + yt .* sin(theta);
    yl = yc - yt .* cos(theta);
end

function [x_c, y_c, x_u, y_u, x_l, y_l, theta] = plot_airfoil(t, m, p)
    x_c = linspace(0, 1, 100);  % chord
    y_c = zeros(size(x_c));
    y_t = zeros(size(x_c));
    x_u = zeros(size(x_c));
    y_u = zeros(size(x_c));
    x_l = zeros(size(x_c));
    y_l = zeros(size(x_c));
    theta = zeros(size(x_c));

    for i = 2:length(x_c)
        y_c(i) = camber_line(m, p, x_c(i));
        y_t(i) = half_thickness(t, x_c(i));
        theta(i) = atan((y_c(i) - y_c(i-1)) / (x_c(i) - x_c(i-1)));
        [x_u(i), y_u(i)] = upper_line(x_c(i), y_t(i), y_c(i), theta(i));
        [x_l(i), y_l(i)] = lower_line(x_c(i), y_t(i), y_c(i), theta(i));
    end
end

function [dy_dxc_theta, A0, A1, c_lift] = thin_airfoil_deta_dx( p, alpha,m)
    % Compute x/c transformation
    theta = linspace(0,pi,20);
    x_c = (1 - cos(theta)) / 2;  % Using x = c/2 * (1 - cos(theta))
    dy_dxc_theta = zeros(length(x_c),1);
    for i = 1:length(x_c)
        dy_dxc_theta(i) = deta_dx(x_c(i),p,m);
    end
    A0 = cal_A0(theta, dy_dxc_theta, alpha);
    A1 = cal_An(theta, dy_dxc_theta, 1);
    alpha_lo_val = alpha_lo(theta, dy_dxc_theta);
    c_lift = (2*pi)*(((deg2rad(alpha))- alpha_lo_val));
    %c_lift = 2*pi*(A0 + (A1/2));
end

function [dy_dx] = deta_dx(x_c, p, m)
    if x_c >= 0 && x_c <= p
        dy_dx = (2*m/(p^2)) * (p - (x_c));
    elseif x_c > p && x_c <= 1
        dy_dx = (2*m/((1-p)^2)) * (p - (x_c));
    else
        dy_dx = 0; % Ensure function returns zero if x_c is out of bounds
    end
end

function [A0] = cal_A0(theta, deta_dx_val, alpha)
    % Ensure theta and deta_dx_val are column vectors
    theta = theta(:);
    deta_dx_val = deta_dx_val(:);

    % Compute integral using trapz
    int_A0 = trapz(theta, deta_dx_val);

    % Compute A0
    A0 = deg2rad(alpha) - ((1/pi) * int_A0);
end

function [An] = cal_An(theta, deta_dx_val, n)
    % Ensure theta and deta_dx_val are column vectors
    theta = theta(:);
    deta_dx_val = deta_dx_val(:);
    
    % Compute the integral term (ensure element-wise multiplication)
    int_var = deta_dx_val .* cos(n * theta);  % Ensure element-wise multiplication
    int_An = trapz(theta, int_var);

    % Compute An
    An = (2/pi) * int_An;
end

function [alpha_lo] = alpha_lo(theta, deta_dx_val)
    % Ensure column vectors
    theta = theta(:);
    deta_dx_val = deta_dx_val(:);

    % Compute integral for alpha_lo
    alpha_lo_var = deta_dx_val .* (cos(theta) - 1);
    alpha_lo_int = trapz(theta, alpha_lo_var);

    % Compute alpha_lo
    alpha_lo = -(1/pi) * alpha_lo_int;
end

function [x_cp, dellcP_q3] = thin_airfoil_dellcP(dy_dxc_theta, alpha, U0, N)
    % Define theta range correctly
    theta = linspace(0, pi, 20);
    theta = theta(:);  % Ensure column vector
    dy_dxc_theta = dy_dxc_theta(:);  % Ensure column vector

    % Compute x_cp (center of pressure locations)
    x_cp = (1/2) * (1 - cos(theta));

    % Preallocate for efficiency
    dellcP_q3 = zeros(length(theta), 1);

    % Compute A0 once
    A0 = cal_A0(theta, dy_dxc_theta, alpha);

    % Compute all An values once
    An_vals = zeros(N, 1);
    for j = 1:N
        An_vals(j) = cal_An(theta, dy_dxc_theta, j);
    end

    % Compute circulation and ΔC_p for each theta(i)
    for i = 1:length(theta)
        % Compute summation term sum(An * sin(n*theta))
        sum_An = sum(An_vals .* sin((1:N)' * theta(i)));

        % Prevent singularity at sin(theta) ≈ 0
        sin_theta_i = sin(theta(i));
        if abs(sin_theta_i) < 1e-6
            sin_theta_i = 1e-6; % Avoid division by near-zero values
        end

        % Compute circulation distribution Gamma(theta)
        gamma = 2 * U0 * (A0 * ((1 + cos(theta(i))) / sin_theta_i) + sum_An);

        % Compute ΔC_p (pressure coefficient difference)
        dellcP_q3(i) = 2 * gamma / U0;
    end
end

function [gamma_new, Cl, Cpr, xp] = Panel_method_new(alpha_deg, Gamma, x1, y1, U0)
    alpha = deg2rad(alpha_deg);
    xdat = [x1,y1];[n1,n2]=size(xdat);
    x=xdat(:,1);
    y=xdat(:,2);
    n=n1-1;

    % Definition of panels and tangent vectors
    for j=1:n
    
        % Panel lengths
        plength(j) = sqrt((x(j+1)-x(j)).^2+(y(j+1)-y(j)).^2);
        
        % Control point
        xp(j)=0.5*(x(j+1)+x(j));
        yp(j)=0.5*(y(j+1)+y(j));
        
        % Tangent vectors
        Tx(j)=-(x(j+1)-x(j))./plength(j); 
        Ty(j)=-(y(j+1)-y(j))./plength(j);
        
        % Normal vectors
        Nx(j)=-Ty(j);
        Ny(j)= Tx(j);
    
    end
    
    % A and B matrices
    for i=1:n
        for j=1:n
           if i == j
               Ux(i,i)=0.0; % Self-induction eq. (19) in note
               Uy(i,i)=0.5; % Self-induction eq. (19) in note
           else
               % Induction at control point 'i' in local system 'j'
               % Eqs. (16) and (18) in note
               sx(i,j)=(xp(i)-xp(j))*Tx(j)+(yp(i)-yp(j))*Ty(j);
               sy(i,j)=(xp(i)-xp(j))*Nx(j)+(yp(i)-yp(j))*Ny(j);
               Ux1(i,j)=log( ((sx(i,j)+0.5*plength(j)).^2 + sy(i,j).^2)/((sx(i,j)-0.5*plength(j)).^2 + sy(i,j).^2) )/(4.*pi);
               Uy1(i,j)=( atan( (sx(i,j)+0.5*plength(j))/sy(i,j) ) - atan( (sx(i,j)-0.5*plength(j))/sy(i,j) ) )/(2.*pi);
              % Induction in global system
               Ux2(i,j)=Ux1(i,j)*Tx(j)-Uy1(i,j)*Ty(j);
               Uy2(i,j)=Ux1(i,j)*Ty(j)+Uy1(i,j)*Tx(j);
              % Induction in local system 'i'
               Ux(i,j)=Ux2(i,j)*Tx(i)+Uy2(i,j)*Ty(i);
               Uy(i,j)=Ux2(i,j)*Nx(i)+Uy2(i,j)*Ny(i);
           end
           % Eq. (24) in note by setting sigma=1
           A(i,j)=Uy(i,j); % Influence coefficients in A matrix (normal velocity)
           B(i,j)=Ux(i,j); % Influence coefficients in B matrix (tangential velocity)
        end
    end
    
    % Boundary conditions
    F = -(Nx.*cos(alpha)+Ny.*sin(alpha));
    
    % Solution of system (solution of eq. (26))
    M = A\F';
    
    for i=1:n
        sum1=0.0;
        sum2=0.0;
        for j=1:n
            sum1=sum1+B(i,j)*M(j,1);
            sum2=sum2+A(i,j)*M(j,1);
        end
        Vt1(i)=sum1;
        Vt(i)=sum1+Tx(i)*cos(alpha)+Ty(i)*sin(alpha); %Tangential velocity (eq. 27)
        Vn(i)=sum2+Nx(i)*cos(alpha)+Ny(i)*sin(alpha); %Check of normal velocity (should be equal to zero)
        Cp(i)=1.-Vt(i).^2; %Pressure coefficient without Kutta condition 
    end
    
    %Vortex distribution
    sum=0.0;
    for i=1:n
        sum=sum+(i-1)*(n-i)*plength(i);
    end
    for i=1:n
        vort(i)=(i-1)*(n-i)/sum; %parabolic vortex distribution  (eq. 37)
    end
    
    C=A\B;
    D=B*C;
    Vrt=(A+D)*vort'; %Rotating onset flow (see eq. 40)
    Cl=2.*Gamma/((abs(max(x)-min(x)))); %Lift coefficient (K-J theorem: Cl=2*Gamma/Vc; V=1)
    Urt=Gamma*Vrt'+Vt; %Tangential velocity including circulation (eq. 41)
    for i=1:n
         Cpr(i)=1.-Urt(i).^2; %Pressure coefficient with Kutta condition
    end 
    gamma_new = -(Vt(1)+Vt(end))/(Vrt(1) + Vrt(end));

end

function [x_c_new, y_c_new] = orient_airfoil(x_u, y_u, x_l, y_l)
    [x_u_sorted, sortIdx] = sort(x_u, 'descend'); % Sort x_c_l in descending order
    y_u_sorted = y_u(sortIdx); % Reorder y_c_l based on the sorting index
    x_c_new = [x_u_sorted'; x_l'];
    y_c_new = [y_u_sorted'; y_l'];

    % Remove duplicate points
    [x_c_new, uniqueIdx] = unique(x_c_new, 'stable'); 
    y_c_new = y_c_new(uniqueIdx);
end

function [x_cp_panel, dell_CP_panel] = cP_Panel(x_act_cp_low, act_cp_low, x_act_cp_high, act_cp_high)
    x_cp_panel = linspace(0, 1, 50); % Increased resolution
    x_cp_panel = x_cp_panel(:);
    dell_CP_panel = zeros(length(x_cp_panel), 1);
    
    % Initialize arrays for interpolated Cp values
    cP_low_interp = zeros(length(x_cp_panel), 1);
    cP_high_interp = zeros(length(x_cp_panel), 1);

    for i = 1:length(x_cp_panel)
        % Interpolating Cp values at each x_cp location
        cP_low_interp(i) = interp1(x_act_cp_low, act_cp_low, x_cp_panel(i), 'pchip');
        cP_high_interp(i) = interp1(x_act_cp_high, act_cp_high, x_cp_panel(i), 'pchip');

        % Ensure positive ΔC_p
        dell_CP_panel(i) = abs(cP_high_interp(i) - cP_low_interp(i));
    end
end

function [x_u_new, y_u_new, x_l_new, y_l_new] = sort_x_y_panel(x_u, y_u, x_l, y_l)
    % Sort upper surface in descending order
    [~, sortIdx_u] = sort(x_u, 'descend');
    x_u_sorted = x_u(sortIdx_u);
    y_u_sorted = y_u(sortIdx_u);
    
    % Sort lower surface in descending order
    [~, sortIdx_l] = sort(x_l, 'descend');
    x_l_sorted = x_l(sortIdx_l);
    y_l_sorted = y_l(sortIdx_l);
    
    % Concatenate first point to close the airfoil shape
    x_u_new = [x_u_sorted, x_u_sorted(1)];
    y_u_new = [y_u_sorted, y_u_sorted(1)];
    x_l_new = [x_l_sorted, x_l_sorted(1)];
    y_l_new = [y_l_sorted, y_l_sorted(1)];
end

function [x_u_sorted, y_u_sorted, CP_u_sorted, x_l_sorted, y_l_sorted, CP_l_sorted] = sort_xfoil_CP(x_u, y_u, CP_u, x_l, y_l, CP_l)
    % Sort upper surface in descending order
    [~, sortIdx_u] = sort(x_u, 'descend');
    x_u_sorted = x_u(sortIdx_u);
    y_u_sorted = y_u(sortIdx_u);
    CP_u_sorted = CP_u(sortIdx_u);
    
    % Sort lower surface in descending order
    [~, sortIdx_l] = sort(x_l, 'descend');
    x_l_sorted = x_l(sortIdx_l);
    y_l_sorted = y_l(sortIdx_l);
    CP_l_sorted = CP_l(sortIdx_l);

end
