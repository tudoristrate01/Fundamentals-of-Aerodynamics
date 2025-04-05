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
[x_c_2312, y_c_2312, x_u_2312, y_u_2312, x_l_2312, y_l_2312, theta_2312, X_Q2_2312, dx_2312] = plot_airfoil(t_2312, m_2312, p_2312);
[x_c_2324, y_c_2324, x_u_2324, y_u_2324, x_l_2324, y_l_2324, theta_2324, X_Q2_2324, dx_2324] = plot_airfoil(t_2324, m_2324, p_2324);
[x_c_4412, y_c_4412, x_u_4412, y_u_4412, x_l_4412, y_l_4412, theta_4412, X_Q2_4412, dx_4412] = plot_airfoil(t_4412, m_4412, p_4412);
[x_c_4424, y_c_4424, x_u_4424, y_u_4424, x_l_4424, y_l_4424, theta_4424, X_Q2_4424, dx_4424] = plot_airfoil(t_4424, m_4424, p_4424);

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

%% Question 2

% Thin airfoil theory

alpha = -10:0.1:15; % Angle of attack in degrees
Cl_2312 = thin_airfoil_lift(m_2312, p_2312, alpha);
Cl_2324 = thin_airfoil_lift(m_2324, p_2324, alpha);
Cl_4412 = thin_airfoil_lift(m_4412, p_4412, alpha);
Cl_4424 = thin_airfoil_lift(m_4424, p_4424, alpha);

figure(3)
plot(alpha, Cl_2312)

% fprintf('Cl for NACA 2312 at %d°: %.3f\n', alpha, Cl_2312);
% fprintf('Cl for NACA 2324 at %d°: %.3f\n', alpha, Cl_2324);
% fprintf('Cl for NACA 4412 at %d°: %.3f\n', alpha, Cl_4412);
% fprintf('Cl for NACA 4424 at %d°: %.3f\n', alpha, Cl_4424);

% Panel Meth


%% Functions

function [yc] = camber_line(m, p, x_c)
    if x_c <= p
        yc = m / p^2 * (2 * p * x_c - x_c.^2);
    else
        yc = m / (1 - p)^2 * (1 - 2*p + 2 * p * x_c - x_c.^2);
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

function [x_c, y_c, x_u, y_u, x_l, y_l, theta, x_Q2, dx] = plot_airfoil(t, m, p)
    x_c = linspace(0, 1, 100);  % chord
    y_c = zeros(size(x_c));
    y_t = zeros(size(x_c));
    x_u = zeros(size(x_c));
    y_u = zeros(size(x_c));
    x_l = zeros(size(x_c));
    y_l = zeros(size(x_c));
    theta = zeros(size(x_c));

    dx = zeros(size(x_c));
    x_Q2 = zeros(size(x_c));
    
    for i = 2:length(x_c)
        y_c(i) = camber_line(m, p, x_c(i));
        y_t(i) = half_thickness(t, x_c(i));
        theta(i) = atan((y_c(i) - y_c(i-1)) / (x_c(i) - x_c(i-1)));
        [x_u(i), y_u(i)] = upper_line(x_c(i), y_t(i), y_c(i), theta(i));
        [x_l(i), y_l(i)] = lower_line(x_c(i), y_t(i), y_c(i), theta(i));

        x_Q2(i) = x_c(i) / 2 * (1 - cos(theta(i)));
        dx(i) = x_c(i) / 2 * sin(theta(i)) * (theta(i) - theta(i-1));
    end
end

function Cl = thin_airfoil_lift(t, m, p, alpha)

    [x_c, y_c, ~, ~, ~, ~, ~, ~, ~] = plot_airfoil(t, m, p);
    alpha = deg2rad(alpha);
    dY_dx = zeros(size(x_c));

    for i = 1:length(dY_dx)
        if x_c(i) <= p
            dY_dx(i) = 2*m / p^2 * (p - x_c);
        else
            dY_dx(i) = 2*m / (1 - p)^2 * (p - x_c);
        end
    end

    
end