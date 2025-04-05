clc
clear 
close all

cd_cl_fixed = 'Cd_fixed.xlsx';
cd_cl_free = 'Cd_free.xlsx';

% Get the sheet names
[~, sheets] = xlsfinfo(cd_cl_fixed);
[~, sheets_free] = xlsfinfo(cd_cl_free);

% Create a figure for the subplots
figure;

% Loop through each sheet and read the data
for i = 1:length(sheets)
    sheetName = sheets{i};
    freeSheetName = sheets_free{i};
    data = readtable(cd_cl_fixed, 'Sheet', sheetName, 'Range', 'A10:G34');
    data_free = readtable(cd_cl_free, 'Sheet', freeSheetName, 'Range', 'A10:G34');
    
    % Create a subplot for each sheet
    subplot(2, 2, i);
    plot(data{:, 3}, data{:, 2}, 'o-');
    hold on
    plot(data_free{:, 3}, data_free{:, 2}, 'o-')
    grid on;
    title(['Plot for NACA', sheetName], 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Drag Coefficient $[-]$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Lift Coefficient $[-]$', 'Interpreter', 'latex', 'FontSize', 16);
    % Customize grid and axis
    ax = gca;
    ax.XMinorGrid = 'on'; % Enable minor grid for better readability
    ax.YMinorGrid = 'on';
    ax.FontSize = 16; % Increase font size for better visibility
    ax.LineWidth = 1.5; % Thicker axis lines

    % Add X and Y Axes at (0,0) but EXCLUDE them from the legend
    xline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % X-axis at y = 0
    yline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Y-axis at x = 0
    
end

% Adjust the figure's PaperPosition property
set(gcf, 'PaperPositionMode', 'auto');

% Save the figure as a PDF using the '-bestfit' option
print(gcf, 'Cl_vs_Cd_Q5.pdf', '-dpdf', '-bestfit');