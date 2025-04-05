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
a = 1.; % Radius of cylinder
U = 1.; % Onset velocity

% Mesh
x = meshgrid(-5:.1:5);
y = x';
% Polar coordinates
r = sqrt(x.^2+y.^2);
theta=atan2(y,x);

% Exclusion of inside-circle points
for n = 1:length(x)
    for m = 1:length(y)
        if r(n,m) <=  a - 5e-3
            r(n,m) = NaN;
            x(n,m) = NaN;
            y(n,m) = NaN;
        end
    end
end

% Aerodynamic potential and stream function
phi = U*r.*cos(theta).*(1+a^2./r.^2);
psi = U*r.*sin(theta).*(1-a^2./r.^2);

% Velocity field
ur = U*cos(theta).*(1-a^2./r.^2);
ut = -U*sin(theta).*(1+a^2./r.^2);
ux = ur.*cos(theta)-ut.*sin(theta);
uy = ur.*sin(theta)+ut.*cos(theta);


% Velocity field around circle
angle = 0.05:.05:2*pi;
x_c = a*cos(angle);
y_c = a*sin(angle);
ut_c = -2*U*sin(angle);
ur_c = 0.;

figure
plot(x_c,-1.+ut_c.^2,'-k','LineWidth',2);
grid
title('Pressure coefficient around cylinder: -Cp');
hold off

figure
plot(angle,-1.+ut_c.^2,'-k','LineWidth',2);
grid
title('Pressure coefficient around cylinder: -Cp');
hold off

figure
hold on
contour(x,y,psi,-10:.1:10)
fill(x_c,y_c,'y')
axis equal
axis([-5 5 -3 3])
title('Streamlines around Circle');
hold off

figure
scale=4.;
quiver(x,y,ux,uy,scale)
hold on
plot(x_c,y_c,'k')
axis equal
axis([-3 3 -3 3])
title('Velocity vectors around Circle');
hold off