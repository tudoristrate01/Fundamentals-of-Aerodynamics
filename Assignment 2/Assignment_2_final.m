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

fig1 = figure;
set(fig1, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
hold on
yline(0, 'k', 'LineWidth', 2);
for i = 1:length(Aspect_Ratio)
    plot (AoA, C_l(i, :), LineWidth=3)
end
hold off
grid on
legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex')
xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Lift coefficient [-]', 'FontSize', 14, 'Interpreter', 'latex')
title('Lift coefficient vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
exportgraphics(fig1, 'Q1_Induced_CL_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');

fig2 = figure;
set(fig2, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
figure (2)
hold on
yline(0, 'k', 'LineWidth', 2);
for i = 1:length(Aspect_Ratio)
    plot (AoA, Cd_i(i, :), LineWidth=3)
end
hold off
grid on
legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Induced Drag coefficent [-]', 'FontSize', 14, 'Interpreter', 'latex')
title('Induced Drag coefficent vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
exportgraphics(fig2, 'Q1_Induced_Cd_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');

fig3 = figure;
set(fig3, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
hold on
yline(0, 'k', 'LineWidth', 2);
for i = 1:length(Aspect_Ratio)
    plot (AoA, CD(i, :), LineWidth=3)
end
hold off
grid on
legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Drag coefficent [-]', 'FontSize', 14, 'Interpreter', 'latex')
title('Drag coefficent vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
exportgraphics(fig3, 'Q1_Cd_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');

fig4 = figure;
set(fig4, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
hold on
yline(0, 'k', 'LineWidth', 2);
for i = 1:length(Aspect_Ratio)
    plot (AoA, alpha_i(i, :), LineWidth=3)
end
hold off
grid on
legend('AR 4', 'AR 6', 'AR 8', 'AR 10', 'AR inf', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
xlabel('Angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Induced angle of attack [Deg]', 'FontSize', 14, 'Interpreter', 'latex')
title('Induced angle of attack vs Angle of attack', 'FontSize', 16, 'Interpreter', 'latex')
exportgraphics(fig4, 'Q1_Induced_AoA_vs_AoA.pdf', 'Resolution', 300, 'ContentType', 'vector');




%% Question 2 Calculations

AoA_Q2 = [0, 5, 10];

m0 = 6.5958;
C = 1;

n = 50;
B = Aspect_Ratio .* C;
theta = linspace(0.01, pi, 180);
x_norm_q2 = linspace(-1,1,180);
A_interm = zeros(n, length(AoA_Q2), length(B));


for i = 1:length(AoA_Q2)
    for j = 1:length(B)
        A_interm(:, i, j) = compute_An_q2(n, B(j), m0, C, AoA_Q2(i), alpha_low, theta);
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

        end
    end
end
x_pos = zeros(length(theta),length(B));
for i = 1:length(B)
    x_pos(:,i) = cos(theta)*(B(i)/2);
end


for i = 1:length(AoA_Q2)
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);
    hold on;

    hLine = yline(0, 'k', 'LineWidth', 2);  % thick horizontal line at y = 0
    hLine.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
    vLine = xline(0, 'k', 'LineWidth', 2);  % thick horizontal line at y = 0
    vLine.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend

    for j = 1:length(B)
        plot(x_norm_q2, rad2deg(squeeze(alpha_i_q2(i, j, :))), 'DisplayName', ['AR = ' num2str(Aspect_Ratio(j))], LineWidth=2);
    end

    title(['Induced AoA vs. Theta for AoA = ' num2str(AoA_Q2(i))], 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('Normalized Spanwise x_{norm} (m)', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Induced AoA (deg)', 'FontSize', 14, 'Interpreter', 'latex');
    legend('Location','best', 'FontSize', 12, 'Interpreter', 'latex');
    grid on;
    hold off;
    exportgraphics(fig, ['Q2_Induced_AoA_vs_Theta' num2str(AoA_Q2(i)) '.pdf'], 'Resolution', 500, 'ContentType', 'vector');
end

%% For lift 

AoA_lift = linspace(-4, 10, 15);
A_1_lift = zeros(1, length(AoA_lift), length(B));
CL_Q2 = zeros(length(AoA_lift), length(B));



for i = 1:length(AoA_lift)
    for j = 1:length(B)
        A_1_lift(:, i, j) = compute_An_q2(1, B(j), m0, C, AoA_lift(i), alpha_low, theta);
        CL_Q2(i,j) = pi*Aspect_Ratio(j)*A_1_lift(:,i,j);
    end
end


fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;
for j = 1:length(B)
    plot(AoA_lift, CL_Q2(:, j), 'DisplayName', ['AR = ' num2str(Aspect_Ratio(j))], 'LineWidth', 2);
end

xlabel('Angle of Attack (deg)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Lift Coefficient (-)', 'FontSize', 14, 'Interpreter', 'latex');
title('C_L vs. AoA for Different AR', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q2_CL_vs_AoA.pdf', 'Resolution', 500, 'ContentType', 'vector');

%% For Drag

n_drag = 50;
A_n_drag = zeros(n_drag, length(AoA_lift), length(B));
dell_Q2 = zeros(length(AoA_lift), length(B));
CD_Q2 = zeros(length(AoA_lift), length(B));

for i = 1:length(AoA_lift)
    for j = 1:length(B)
        A_n_drag(:, i, j) = compute_An_q2(n_drag, B(j), m0, C, AoA_lift(i), alpha_low, theta);
        
        dell_sum = 0;  % Reset dell_sum for each iteration of i or j
        for k = 2:n_drag
            dell_const = k * ((A_n_drag(k,i,j) / A_n_drag(1,i,j))^2);
            dell_sum = dell_sum + dell_const;
        end
        
        dell_Q2(i,j) = dell_sum;
        CD_Q2(i,j) = ((CL_Q2(i,j)^2) / (pi * Aspect_Ratio(j))) * (1 + dell_Q2(i,j));
    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;
for j = 1:length(B)
    plot(AoA_lift, CD_Q2(:, j), 'DisplayName', ['AR = ' num2str(Aspect_Ratio(j))], 'LineWidth', 2);
end

xlabel('Angle of Attack (deg)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Drag Coefficient (-)', 'FontSize', 14, 'Interpreter', 'latex');
title('C_D vs. AoA for Different AR', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q2_CD_vs_AoA.pdf', 'Resolution', 500, 'ContentType', 'vector');


%% Question 3

TR = [0.2, 0.44, 0.6, 0.8, 1.0];
AR_3 = 6;
AoA_Q3 = 6; % In deg
norm_x = linspace(-1,1,500);
% Considering c_root as 1
c_q3 = zeros(length(norm_x), length(TR));
b_q3 = zeros(length(TR),1);
c_cap = zeros(length(TR),1);
x_q3 = zeros(length(norm_x), length(TR));
theta_q3 = zeros(length(norm_x), length(TR));
n_q3 = 100;
for i = 1:length(norm_x)
    for j = 1:length(TR)
        c_q3(i,j) = 1-((1-TR(j))*abs(norm_x(i)));
        b_q3(j) = (AR_3*1*(1+TR(j)))/2;
        x_q3(i,j) = norm_x(i)*(b_q3(j)/2);
        theta_q3(i,j) = acos(-2*x_q3(i,j)/b_q3(j));  % Note the negative sign
        c_cap(j) = 1*((1+TR(j))/2);
    end
end

theta_q3_deg = rad2deg(theta_q3);
A_q3 = zeros(n_q3, length(TR));




for j = 1:length(TR)
    A_q3(:, j) = compute_An_q3(n_q3, b_q3(j), m0, c_q3(:,j), AoA_Q3, alpha_low, theta_q3(:,j));
end


for i = 1:length(theta_q3(:, 1))  % Loop over rows of theta_q3
    for k = 1:length(TR)
        Gamma_sum_q3 = 0;
        for j = 1:n_q3
            Gamma_inter_q3(j, i) = A_q3(j, k) * sin(j * theta_q3(i, k));  % Fixed indexing
            Gamma_sum_q3 = Gamma_inter_q3(j, i) + Gamma_sum_q3;
        end
        
        Gamma_q3(i, k) = Gamma_sum_q3 * 2 * b_q3(k) / c_cap(k);  % Corrected indexing
    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on; 
for k = 1:length(TR)
    plot(norm_x(30:470), Gamma_q3(30:470,k) , 'DisplayName', ['TR = ' num2str(TR(k))], 'LineWidth', 2);
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Circulation', 'FontSize', 14, 'Interpreter', 'latex');
title('Circulation vs. Normalized Spanwise Position for Different TR', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q3_Circulation_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');



fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;
for k = 1:length(TR)
    plot(norm_x(30:470), c_q3(30:470,k) , 'DisplayName', ['TR = ' num2str(TR(k))], 'LineWidth', 2);
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Chord length', 'FontSize', 14, 'Interpreter', 'latex');
title('Chord vs. Normalized Spanwise Position for Different TR', 'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q3_Chord_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');



%% Induced AoA

alpha_i_q3 = zeros(length(TR), length(theta_q3));
for k = 1:length(theta_q3)-1
    for j = 1:length(B)
        alpha_interm_sum_q3 = 0;  % Reset sum accumulator for each (i,j,k) combination
        for z = 1:n_q3
            alpha_interm_q3 = z * A_q3(z,j) * (sin(z*theta_q3(k)) / sin(theta_q3(k)));
            alpha_interm_sum_q3 = alpha_interm_sum_q3 + alpha_interm_q3;  % Add to the sum
        end
        alpha_i_q3(j,k) = alpha_interm_sum_q3;  % Assign the final sum after the loop over z

    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

% hLine = yline(0, 'k', 'LineWidth', 2);
% hLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
vLine = xline(0, 'k', 'LineWidth', 2);
vLine.Annotation.LegendInformation.IconDisplayStyle = 'off';

hold on;
for j = 1:length(TR)
    plot(norm_x(30:470), rad2deg(alpha_i_q3(j, 30:470)), 'DisplayName', ['TR = ' num2str(TR(j))], 'LineWidth', 2);
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Induced Angle of Attack (deg)', 'FontSize', 14, 'Interpreter', 'latex');
title('Induced Angle of Attack vs. Normalized Spanwise Position for Different TR', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q3_Ind_AoA_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');

%% Local lift coefficient
cl_q3 = zeros(length(TR), length(theta_q3));
for k = 1:length(theta_q3)
    for j = 1:length(b_q3)
        cl_interm_sum_q3 = 0;  % Reset sum accumulator for each (i,j,k) combination
        for z = 1:n_q3
            cl_interm_q3 = A_q3(z,j) * (sin(z*theta_q3(k)));
            cl_interm_sum_q3 = cl_interm_sum_q3 + cl_interm_q3;  % Add to the sum
        end
        cl_q3(j,k) = (4*b_q3(j)/c_q3(k,j))*cl_interm_sum_q3;  % Assign the final sum after the loop over z

    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

% hLine = yline(0, 'k', 'LineWidth', 2);
% hLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
vLine = xline(0, 'k', 'LineWidth', 2);
vLine.Annotation.LegendInformation.IconDisplayStyle = 'off';

hold on;
for j = 1:length(TR)
    plot(norm_x(30:470), cl_q3(j, 30:470), 'DisplayName', ['TR = ' num2str(TR(j))], 'LineWidth', 2);  % Plot alpha_i_q3 for each TR
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Local lift coefficient (-)', 'FontSize', 14, 'Interpreter', 'latex');
title('Local lift coefficient vs. Normalized Spanwise Position for Different TR', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q3_Local_CL_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');


%% Local drag coefficient

cd_q3 = zeros(length(TR), length(theta_q3));
for k = 1:length(theta_q3)
    for j = 1:length(B)
        cd_q3(j,k) = cl_q3(j,k)*alpha_i_q3(j,k);
    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

% hLine = yline(0, 'k', 'LineWidth', 2);
% hLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
vLine = xline(0, 'k', 'LineWidth', 2);
vLine.Annotation.LegendInformation.IconDisplayStyle = 'off';

hold on;
for j = 1:length(TR)
    plot(norm_x(30:470), cd_q3(j, 30:470), 'DisplayName', ['TR = ' num2str(TR(j))], 'LineWidth', 2);  % Plot alpha_i_q3 for each TR
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Local drag coefficient (-)', 'FontSize', 14, 'Interpreter', 'latex');
title('Local drag coefficient vs. Normalized Spanwise Position for Different TR', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q3_Local_CD_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');





%% Question 4

AR_4 = 6;
% Considering C = 1
B_4 = AR_4*1;
alpha_g_tip = [0; 2; 4; 6; 8];
alpha_g0 = 4;
theta_q4 = linspace(0.01,pi,500);
x_q4 = (cos(theta_q4)*(B_4/2));
x_q4_plot = x_q4/(B_4/2);
alpha_geo_q4 = zeros(length(x_q4), length(alpha_g_tip));
n_4 = 200;


for j = 1:length(x_q4)
    for i = 1:length(alpha_g_tip)
        alpha_geo_q4(j,i) = deg2rad(alpha_g0) + ((deg2rad(alpha_g_tip(i))-deg2rad(alpha_g0))*abs(x_q4(j)/(B_4/2)));
    end
end
A_q4 = zeros(n_4, length(alpha_g_tip));
for j = 1:length(alpha_g_tip)
    A_q4(:, j) = compute_An_q4(n_4, B_4, m0, 1, rad2deg(alpha_geo_q4(:,j)), alpha_low, theta_q4);
end


for i = 1:length(theta_q4)  % Loop over rows of theta_q3
    for k = 1:length(alpha_g_tip)
        Gamma_sum_q4 = 0;
        for j = 1:n_4
            Gamma_inter_q4 = A_q4(j, k) * sin(j * theta_q4(i));  % Fixed indexing
            Gamma_sum_q4 = Gamma_inter_q4 + Gamma_sum_q4;
        end
        
        Gamma_q4(i, k) = Gamma_sum_q4 * 2 * B_4 / 1;  % Corrected indexing
    end
end


fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;

for k = 1:length(alpha_g_tip)
    plot(x_q4_plot(20:480), alpha_geo_q4(20:480, k), 'DisplayName', sprintf('$\\alpha_{\\mathrm{tip}} = %.2f$', alpha_g_tip(k)), 'LineWidth', 2);  
end

xlabel('Normalized Spanwise Position)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Alpha geo', 'FontSize', 14, 'Interpreter', 'latex');
title('Geometric angle vs. Normalized Spanwise Position for Different angles', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q4_Alpha_geo_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');



fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;

for k = 1:length(alpha_g_tip)
    plot(x_q4_plot(20:480), Gamma_q4(20:480, k), 'DisplayName', sprintf('$\\alpha_{\\mathrm{tip}} = %.2f$', alpha_g_tip(k)), 'LineWidth', 2);  
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Circulation', 'FontSize', 14, 'Interpreter', 'latex');
title('Circulation vs. Normalized Spanwise Position for Different angles', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q4_Circulation_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');

%% Induced AoA 

% Ensure correct dimensions
alpha_i_q4 = zeros(length(alpha_g_tip), length(theta_q4));



for j = 1:length(alpha_g_tip)
    for i = 1:length(theta_q4)
        alpha_q4_sum = 0;
        for k = 1:n_4
            alpha_q4_sum = alpha_q4_sum + (k * A_q4(k, j) * sin(k * theta_q4(i))) / sin(theta_q4(i));
        end
        alpha_i_q4(j, i) = alpha_q4_sum;  % Store final sum
    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;

for j = 1:length(alpha_g_tip)
    plot(x_q4_plot(20:480), (rad2deg(alpha_i_q4(j, 20:480))), 'DisplayName', sprintf('$\\alpha_{\\mathrm{tip}} = %.2f$', alpha_g_tip(j)), 'LineWidth', 2);
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Induced Angle of Attack (deg)', 'FontSize', 14, 'Interpreter', 'latex');
title('Induced Angle of Attack vs. Normalized Spanwise Position for Different angles', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q4_induced_AoA_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');



%% Local lift coefficient
cl_q4 = zeros(length(alpha_g_tip), length(theta_q4));
for k = 1:length(theta_q4)
    for j = 1:length(alpha_g_tip)
        cl_interm_sum_q4 = 0;  % Reset sum accumulator for each (i,j,k) combination
        for z = 1:n_4
            cl_interm_q4 = A_q4(z,j) * (sin(z*theta_q4(k)));
            cl_interm_sum_q4 = cl_interm_sum_q4 + cl_interm_q4;  % Add to the sum
        end
        cl_q4(j,k) = (4*B_4/1)*cl_interm_sum_q4;  % Assign the final sum after the loop over z

    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;
for j = 1:length(TR)
    plot(x_q4_plot(20:480), cl_q4(j, 20:480), 'DisplayName', sprintf('$\\alpha_{\\mathrm{tip}} = %.2f$', alpha_g_tip(j)), 'LineWidth', 2);
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Local lift coefficient (-)', 'FontSize', 14, 'Interpreter', 'latex');
title('Local lift coefficient vs. Normalized Spanwise Position for Different angles', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q4_local_CL_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');


%% Local drag coefficient

cd_q4 = zeros(length(alpha_g_tip), length(theta_q4));
for k = 1:length(theta_q4)
    for j = 1:length(alpha_g_tip)
        cd_q4(j,k) = cl_q4(j,k)*alpha_i_q4(j,k);
    end
end

fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100, 100, 800, 600]);

hold on;  % Hold the plot to overlay multiple plots for different TR values
for j = 1:length(alpha_g_tip)
    plot(x_q4_plot(20:480), cd_q4(j, 20:480), 'DisplayName', sprintf('$\\alpha_{\\mathrm{tip}} = %.2f$', alpha_g_tip(j)), 'LineWidth', 2);
end

xlabel('Normalized Spanwise Position', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Local drag coefficient (-)', 'FontSize', 14, 'Interpreter', 'latex');
title('Local drag coefficient vs. Normalized Spanwise Position for Different angles', 'FontSize', 16, 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
hold off;
exportgraphics(fig, 'Q4_Llocal_CD_vs_Norm_x.pdf', 'Resolution', 500, 'ContentType', 'vector');




%% Functions

function [C_l, Cd_i, alpha_i_rad] = calculate_Q1(AR, alpha, alpha_low)  % eliptic wing formula
    C_l = 2*pi / (1 + 2 / AR) * deg2rad(alpha - alpha_low);
    Cd_i = 1 / (pi * AR) * C_l^2;
    alpha_i_rad = deg2rad(alpha - alpha_low) / (1 + AR/2);
end

function An = compute_An_q2(n, b, m0, C, alpha, alpha_L0, theta)
    % SOLVE_COEFFICIENTS Solves for the first n coefficients A_n
    % Inputs:
    %   n        - Number of coefficients to compute
    %   b        - Given constant parameter
    %   m0       - Function handle or vector of m0(theta)
    %   C        - Array of values corresponding to chord length at each theta
    %   alpha    - Function handle or vector of alpha(theta)
    %   alpha_L0 - Function handle or vector of alpha_L0(theta)
    %   theta    - Vector of theta values for numerical evaluation

    % Construct the system of equations
    S = zeros(length(theta), n);
    T = zeros(length(theta), n);
    
    for k = 1:n
        for j = 1:length(theta)
            if theta(j) == 0
                theta(j) = 1e-6;
            end
            % if theta(j) == pi
            %     theta(j) = 0.01;
            % end
            S(j, k) = (-4*b*sin(k*theta(j))) / (m0 * C);
            T(j, k) = (k*sin(k*theta(j))/sin(theta(j)));
        end
    end
    
    % Right-hand side of the equation
    RHS = zeros(length(theta),1);
    for i = 1:length(theta)
        RHS(i) = -deg2rad(alpha - alpha_L0);
    end

    
    % Solve the system
    A_matrix = S - T;
    An = A_matrix \  RHS;

end

function An = compute_An_q3(n, b, m0, C, alpha, alpha_L0, theta)
     % SOLVE_COEFFICIENTS Solves for the first n coefficients A_n
    % Inputs:
    %   n        - Number of coefficients to compute
    %   b        - Given constant parameter
    %   m0       - Function handle or vector of m0(theta)
    %   C        - Function handle or vector of C(theta)
    %   alpha    - Function handle or vector of alpha(theta)
    %   alpha_L0 - Function handle or vector of alpha_L0(theta)
    %   theta    - Vector of theta values for numerical evaluation

    
    % Construct the system of equations
    S = zeros(length(theta), n);
    T = zeros(length(theta), n);
    
    for k = 1:n
        for j = 1:length(theta)
            if theta(j) == 0
                theta(j) = 1e-6;
            end
            % if theta(j) == pi
            %     theta(j) = 0.01;
            % end
            S(j, k) = (-4*b*sin(k*theta(j))) / (m0 * C(j));
            T(j, k) = (k*sin(k*theta(j))/sin(theta(j)));
        end
    end
    
    % Right-hand side of the equation
    RHS = zeros(length(theta),1);
    for i = 1:length(theta)
        RHS(i) = -deg2rad(alpha - alpha_L0);
    end

    
    % Solve the system
    A_matrix = S - T;
    An = A_matrix \  RHS;

end

function An = compute_An_q4(n, b, m0, C, alpha, alpha_L0, theta)
    % SOLVE_COEFFICIENTS Solves for the first n coefficients A_n
    % Inputs:
    %   n        - Number of coefficients to compute
    %   b        - Given constant parameter
    %   m0       - Function handle or vector of m0(theta)
    %   C        - Function handle or vector of C(theta)
    %   alpha    - Function handle or vector of alpha(theta)
    %   alpha_L0 - Function handle or vector of alpha_L0(theta)
    %   theta    - Vector of theta values for numerical evaluation

    
    % Construct the system of equations
    S = zeros(length(theta), n);
    T = zeros(length(theta), n);
    
    for k = 1:n
        for j = 1:length(theta)
            if theta(j) == 0
                theta(j) = 0.01;
            end
            % if theta(j) == pi
            %     theta(j) = 0.01;
            % end
            S(j, k) = (-4*b*sin(k*theta(j))) / (m0 * C);
            T(j, k) = (k*sin(k*theta(j))/sin(theta(j)));
        end
    end
    
    % Right-hand side of the equation
    RHS = zeros(length(theta),1);
    for i = 1:length(theta)
        RHS(i) = -deg2rad(alpha(i) - alpha_L0);
    end

    
    % Solve the system
    A_matrix = S - T;
    An = A_matrix \  RHS;

end

