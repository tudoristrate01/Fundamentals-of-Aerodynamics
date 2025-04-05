clc
close all
clear all

values = load('table_values.txt');

alpha = values(:,1);
c_l = values(:,2);
c_d = values(:,3);
c_m_c4 = values(:,4);

x_div_c = 1/4 - c_m_c4 ./ c_l;

% figure(1)
% plot(c_d, c_l)
% grid on
% xlabel('cd')
% ylabel('cl')

figure(2)
% loglog(alpha, x_div_c)
plot(alpha, x_div_c)
grid on
xlabel('alpha')
ylabel('x_c')