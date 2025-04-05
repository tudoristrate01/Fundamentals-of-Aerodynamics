%-------------------------------------------------------------------------
%                  Plot of flow around cylinder
%
%                     Jens N. Sørensen
%                DTU Wind Energy - 05-02-2024
%-------------------------------------------------------------------------

clear all
close all
clc

% Input data
a = 3.; % Radius of cylinder
U = 1.; % Onset velocity
Gamma = 20.; % Circulation around the cylinder (adjust for different rotation rates)

% Mesh
[x, y] = meshgrid(-20:.1:20);  % change here the margins of the figures
r = sqrt(x.^2 + y.^2);
theta = atan2(y, x);

% Exclusion of inside-circle points
r(r <= a - 5e-3) = NaN;
x(isnan(r)) = NaN;
y(isnan(r)) = NaN;

% Aerodynamic potential and stream function with circulation
phi = U * r .* cos(theta) .* (1 + a^2 ./ r.^2) + (Gamma / (2 * pi)) * theta;
psi = U * r .* sin(theta) .* (1 - a^2 ./ r.^2) - (Gamma / (2 * pi)) * log(r);

% Velocity field with circulation
ur = U * cos(theta) .* (1 - a^2 ./ r.^2);
ut = -U * sin(theta) .* (1 + a^2 ./ r.^2) + Gamma ./ (2 * pi * r);
ux = ur .* cos(theta) - ut .* sin(theta);
uy = ur .* sin(theta) + ut .* cos(theta);

% Velocity field around cylinder (on the surface)
angle = 0.05:.05:2*pi;
x_c = a * cos(angle);
y_c = a * sin(angle);
ut_c = -2 * U * sin(angle) + Gamma / (2 * pi * a);
ur_c = 0; % No radial component on the surface

% Pressure coefficient around cylinder: Cp = 1 - (V/U)^2
Cp = 1 - (ut_c / U).^2;

% Plot Pressure Coefficient vs angle
figure
plot(angle, Cp, '-k', 'LineWidth', 2);
grid on
title('Pressure Coefficient around Rotating Cylinder');
xlabel('\theta (rad)');
ylabel('C_p');
axis([0 2*pi -4 2]);
hold off

% Plot Pressure Coefficient along the cylinder surface
figure
plot(x_c, Cp, '-k', 'LineWidth', 2);
grid on
title('Pressure Coefficient Distribution on Cylinder Surface');
xlabel('x');
ylabel('C_p');
axis equal
hold off

% Plot Streamlines
figure
hold on
contour(x, y, psi, -10:.1:10)
fill(x_c, y_c, 'y') % Cylinder filled for visualization
axis equal
axis([-5 5 -3 3])
title('Streamlines around Rotating Cylinder');
hold off

% Plot Velocity Vectors
figure
scale = 4.;
quiver(x, y, ux, uy, scale)
hold on
plot(x_c, y_c, 'k', 'LineWidth', 2)
axis equal
axis([-3 3 -3 3])
title('Velocity Vectors around Rotating Cylinder');
hold off