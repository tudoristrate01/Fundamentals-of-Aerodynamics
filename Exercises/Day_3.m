clear all 
close all 
clc

% Input data 
alphad = 30.;  % Angle of attack 
alpha = alphad*pi/180.;
a = 1.0;  % Radius of cylinder (a>c) 
ros = 0.25;  % ratio of semi-axes 
c = a*sqrt((1-ros)/(1+ros));  % Transformationconstant (zeta=z+c^2/z) 
c2 = c^2;

major_str = num2str(a+c2/a);  % Major axis 
minor_str = num2str(a-c2/a);  % Minor axis 
theta_str = num2str(alphad);

% Mesh 
x = meshgrid(-5:.1:5);
y = x';

% Complex mesh plane 
z = x + 1i*y;

% Exclusionof inside-circle points 
for n = 1:length(x) 
    for m = 1:length(y) 
        if abs(z(n,m)) <=  a - 5e-4 
        z(n,m) = NaN;
        end 
    end 
end 

% Aerodynamic potential 
f = exp(-1i*alpha).*z + (exp(1i*alpha)*a^2)./z;

% Joukowsky Transformation
zeta = z + c2./z;

% Transformed potential;Is actually not needed or used 
fz=zeta*exp(-1i*alpha)+((a/c)^2*exp(1i*alpha)-exp(-1i*alpha))*(zeta/2-sqrt((zeta/2).^2-c^2));

% Velocity field around circle 
dfdz = exp(-1i*alpha) - (exp(1i*alpha)*a^2)./z.^2;
u1 = real(dfdz);
v1 = -imag(dfdz);

% Velocity field around ellipse 
W = dfdz./(1.-c2./z.^2);
u = real(W);
v = - imag(W);

% Contour of Circle and Ellipse 
angle = 0.0:.05:2*pi+.05;
z_circle = a*(cos(angle)+1i*sin(angle));
z_ellipse = z_circle+c2./z_circle;

% Check of transformationoncontour around the ellipse 
fcc = exp(-1i*alpha).*z_circle + (exp(1i*alpha)*a^2)./z_circle;
fee=z_ellipse*exp(-1i*alpha)+((a/c)^2*exp(1i*alpha)-exp(-1i*alpha))*(z_ellipse/2 .* sqrt((z_ellipse/2).^2-c^2));
dff=fcc-fee;  % Error function(should be zero to machine accuracy)

% Circle: Velocity oncontour 
dfcdz = exp(-1i*alpha) - (exp(1i*alpha)*a^2)./z_circle.^2;
uc1 = real(dfcdz);
vc1 = -imag(dfcdz);
uctheta = -uc1.*sin(angle)+vc1.*cos(angle);
wc = sqrt(uc1.^2+vc1.^2);

% Ellipse: Velocity oncontour 
Wa = dfcdz./(1.-c2./z_circle.^2);
ua = real(Wa);
va = - imag(Wa);
wa = sqrt(ua.^2+va.^2);

% Exploit sign of circle flow to define sign of ellipse flow oncontour 
for i=1:length(angle)
    watheta(i)=sign(uctheta(i))*wa(i);
end

% Plots of solutions 
figure (1)
plot(real(z_ellipse),imag(z_ellipse),'-b','LineWidth',2,'MarkerSize',8);
hold on
plot(real(z_circle),imag(z_circle),'-r','LineWidth',2,'MarkerSize',8);
axis equal 
title(strcat('Ellipse.  Major axis:',major_str,';Minor axis:',minor_str));
grid 
hold off

figure (2)
plot(angle,watheta,'-b',angle,uctheta,'-r','LineWidth',2);
grid 
xlabel('Angle (radians)');
ylabel('U(theta) [-]') 
title('Azimuthal velocity distributions');
legend('Ellipse','Cylinder','Location','East');
hold off

figure (3)
plot(angle,1.-(ua.^2+va.^2),'-b',angle,1.-(uc1.^2+vc1.^2),'-r','LineWidth',2);
title(strcat('Pressure coefficients'));
legend('Ellipse','Cylinder','Location','East');
xlabel('Angle (radians)');
ylabel('Cp [-]') 
grid 
hold off

figure (4)
hold on
contour(real(z),imag(z),imag(f),-5:.05:5) 
fill(real(z_circle),imag(z_circle),'y') 
axis equal 
axis([-4 4 -2 2]) 
title('Streamlines around the Circle');
hold off

figure (5)
hold on
contour(real(z),imag(z),real(f),-5:.05:5) 
fill(real(z_circle),imag(z_circle),'y') 
axis equal 
axis([-4 4 -2 2]) 
title('Iso-velocity potential lines around the Circle');
hold off

figure (6)
hold on
contour(real(zeta),imag(zeta),imag(f),-8:.05:8) 
fill(real(z_ellipse),imag(z_ellipse),'y') 
axis equal 
axis([-4 4 -2 2]) 
title('Streamlines around the Ellipse');
hold off

figure (7)
hold on
contour(real(zeta),imag(zeta),real(f),-8:.05:8) 
fill(real(z_ellipse),imag(z_ellipse),'y') 
axis equal 
axis([-4 4 -2 2]) 
title('Iso-velocity potential lines around the Ellipse');
hold off

figure (8)
scale=3.;
quiver(real(z),imag(z),u1,v1,scale) 
hold on
plot(real(z_circle),imag(z_circle),'k') 
axis equal 
axis([-3 3 -2 2]) 
title('Velocity vectors around the Circle');
hold off

figure (9)
scale=3.;
quiver(real(zeta),imag(zeta),u,v,scale) 
hold on
plot(real(z_ellipse),imag(z_ellipse),'k','linewidth',1) 
axis equal 
axis([-3 3 -2 2]) 
title('Velocity vectors around the Ellipse');
hold off