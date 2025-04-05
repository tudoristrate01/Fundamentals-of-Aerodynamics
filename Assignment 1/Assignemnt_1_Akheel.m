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

%% Question 2

% Thin airfoil theory
alpha = linspace(-10, 15, 100);  % Creates 100 points between -10 and 15

% NACA 2312
[x_c_theta_2312, dy_dxc_theta_2312, int_A0_2312, int_A1_2312] = thin_airfoil_deta_dx(x_c_2312, p_2312, m_2312, theta_2312);
[A0_2312, A1_2312, c_lift_2312] = thin_airfoil_clift(alpha, int_A0_2312, int_A1_2312);

% NACA 2324
[x_c_theta_2324, dy_dxc_theta_2324, int_A0_2324, int_A1_2324] = thin_airfoil_deta_dx(x_c_2324, p_2324, m_2324, theta_2324);
[A0_2324, A1_2324, c_lift_2324] = thin_airfoil_clift(alpha, int_A0_2324, int_A1_2324);

% NACA 4412
[x_c_theta_4412, dy_dxc_theta_4412, int_A0_4412, int_A1_4412] = thin_airfoil_deta_dx(x_c_4412, p_4412, m_4412, theta_4412);
[A0_4412, A1_4412, c_lift_4412] = thin_airfoil_clift(alpha, int_A0_4412, int_A1_4412);

% NACA 4424
[x_c_theta_4424, dy_dxc_theta_4424, int_A0_4424, int_A1_4424] = thin_airfoil_deta_dx(x_c_4424, p_4424, m_4424, theta_4424);
[A0_4424, A1_4424, c_lift_4424] = thin_airfoil_clift(alpha, int_A0_4424, int_A1_4424);

figure;
subplot(2,2,1);
plot(alpha, c_lift_2312, 'b-');
xlabel('Angle of attack (deg)');
ylabel('Lift coefficient [-]');
legend('Thin airfoil theory');
grid on
title('Lift coeficient Vs. AoA - NACA 2312');

subplot(2,2,2);
plot(alpha, c_lift_2324, 'b-');
xlabel('Angle of attack (deg)');
ylabel('Lift coefficient [-]');
legend('Thin airfoil theory');
grid on
title('Lift coeficient Vs. AoA - NACA 2324');

subplot(2,2,3);
plot(alpha, c_lift_4412, 'b-');
xlabel('Angle of attack (deg)');
ylabel('Lift coefficient [-]');
legend('Thin airfoil theory');
grid on
title('Lift coeficient Vs. AoA - NACA 4412');

subplot(2,2,4);
plot(alpha, c_lift_4424, 'b-');
xlabel('Angle of attack (deg)');
ylabel('Lift coefficient [-]');
legend('Thin airfoil theory');
grid on
title('Lift coeficient Vs. AoA - NACA 4424');

% Panel Meth

% XFOIL

%% Question 3

% Thin Airfoil theory

rho = 1.225;
U0 = 8;
alpha_q3 = 10;

% NACA 2312
dellcP_q3_2312 = thin_airfoil_dellcP(x_c_2312, dy_dxc_theta_2312, rho, U0, alpha_q3, theta_2312);

% NACA 2324
dellcP_q3_2324 = thin_airfoil_dellcP(x_c_2324, dy_dxc_theta_2324, rho, U0, alpha_q3, theta_2324);

% NACA 4412
dellcP_q3_4412 = thin_airfoil_dellcP(x_c_4412, dy_dxc_theta_4412, rho, U0, alpha_q3, theta_4412);

% NACA 4424
dellcP_q3_4424 = thin_airfoil_dellcP(x_c_4424, dy_dxc_theta_4424, rho, U0, alpha_q3, theta_4424);

figure;

% NACA 2312
subplot(2,2,1);
plot(x_c_2312, dellcP_q3_2312, 'b-');
xlabel('x/c');
ylabel('$\Delta C_P$', 'Interpreter', 'latex');
title('$\Delta C_P$ vs. $x/c$ - NACA 2312', 'Interpreter', 'latex');
grid on
legend('Thin airfoil theory');

% NACA 2324
subplot(2,2,2);
plot(x_c_2324, dellcP_q3_2324, 'b-');
xlabel('x/c');
ylabel('$\Delta C_P$', 'Interpreter', 'latex');
title('$\Delta C_P$ vs. $x/c$ - NACA 2324', 'Interpreter', 'latex');
grid on
legend('Thin airfoil theory');

% NACA 4412
subplot(2,2,3);
plot(x_c_4412, dellcP_q3_4412, 'b-');
xlabel('x/c');
ylabel('$\Delta C_P$', 'Interpreter', 'latex');
title('$\Delta C_P$ vs. $x/c$ - NACA 4412', 'Interpreter', 'latex');
grid on
legend('Thin airfoil theory');

% NACA 4424
subplot(2,2,4);
plot(x_c_4424, dellcP_q3_4424, 'b-');
xlabel('x/c');
ylabel('$\Delta C_P$', 'Interpreter', 'latex');
title('$\Delta C_P$ vs. $x/c$ - NACA 4424', 'Interpreter', 'latex');
grid on
legend('Thin airfoil theory');

function[dellcP_q3] = thin_airfoil_dellcP(x_c, dy_dxc_theta,rho, U0, alpha_q3, theta)
    A0_q3 = zeros(length(x_c),1);
    A1_q3 = zeros(length(x_c),1);
    gamma_q3 = zeros(length(x_c),1);
    pu_pl = zeros(length(x_c),1);
    dellcP_q3 = zeros(length(x_c),1);
    for i = 1:length(x_c)
        A0_q3(i) = deg2rad(alpha_q3) - ((1/pi)*dy_dxc_theta(i));
        A1_q3(i) = (2/pi)*(dy_dxc_theta(i)*cos(theta(i)));
        gamma_q3(i) = U0*pi*(A0_q3(i)+(A1_q3(i)/2));
        pu_pl(i) = -(rho*U0*gamma_q3(i));
        dellcP_q3(i) = pu_pl(i)/(0.5*rho*U0^2);
    end
end



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

function[x_c_theta, dy_dxc_theta, int_A0, int_A1] = thin_airfoil_deta_dx(x_c, p, m, theta)
    x_c_theta = zeros(length(x_c),1);
    dy_dxc_theta = zeros(length(x_c),1);
    dy_dxc_theta_A1 = zeros(length(x_c),1);
    for i = 2:length(x_c)
        x_c_theta(i) = (1/2)*(1-cos(theta(i)));
        if x_c(i)>=0 && x_c(i)<=p
            dy_dxc_theta(i) = ((2*m)/(p^2))*(p - x_c(i));
        elseif x_c(i) > p && x_c(i) <= 1
            dy_dxc_theta(i) = ((2*m)/((1-p)^2))*(p - x_c(i));
        end
        dy_dxc_theta_A1(i) = dy_dxc_theta(i)*cos(theta(i));
    end
    int_A0 = trapz(dy_dxc_theta,theta);
    int_A1 = trapz(dy_dxc_theta_A1,theta);
end

function[A0, A1, c_lift] = thin_airfoil_clift(alpha, int_A0, int_A1)
    A0 = zeros(length(alpha),1);
    A1 = zeros(length(alpha),1);
    c_lift = zeros(length(alpha),1);
    for i = 1:length(alpha)
        A0(i) = deg2rad(alpha(i)) - ((1/pi)*int_A0);
        A1(i) = (2/pi)*int_A1;
        c_lift(i) = 2*pi*(A0(i) + (A1(i)/2));
    end
end