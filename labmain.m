close all
clear
% propriedades do ar (21ºC)
pT = 100658;
rho_air = 1.201;    % Density of air [kg/m^3]
nu_air = 1.515e-5;  % Dynamic viscosity of air
mu_air = 1.819e-5;  % Kinematic viscosity of air [m^2/s]

% dados enunciado
rho_fm = 825;
g = 9.81;
L1 = 1.7;
L2 = 0.155;
L3 = 0.015; % offset da transição
largura = 0.235;
beta = 12.9;
l0 = 0.1;
p0 = rho_fm*g*l0*sind(beta);
%p0 = 325;
latm = 0.460;

lpsref = 0.366;
lpTref = 0.187;
manometro  = [0.360; 0.378; 0.392; 0.406];
l = latm - manometro;
% ponto de referência dentro do ventilador
pTref = (latm - lpTref) * rho_fm * g * sind(beta);
psref = (latm - lpsref) * rho_fm * g * sind(beta);
qref = pTref - psref;
Ueref = sqrt(2 * qref / rho_air);
pslocal = l * rho_fm * g * sind(beta);
qlocal = pTref - pslocal;
Uelocal = sqrt(2 * qlocal / rho_air);
Cp = (pslocal - psref) / qref;

%%
x = [0.225, 0.475, 0.725, 0.975];
save('x', 'x');
Cfabaco = [0.00395, 0.00375, 0.00353, 0.00339]; % Encontrados à mão para critério no calculo dos parametros
save('Cfabaco', 'Cfabaco');

pxraw1 = table2array(readtable('X1_20240430_181354.xlsx', 'Sheet', 'Motor X - Line 1','VariableNamingRule','preserve'));
pxraw2 = table2array(readtable('X2_20240430_175840.xlsx', 'Sheet', 'Motor X - Line 1','VariableNamingRule','preserve'));
pxraw3 = table2array(readtable('X3_20240430_174255.xlsx', 'Sheet', 'Motor X - Line 1','VariableNamingRule','preserve'));
pxraw4 = table2array(readtable('X4_20240430_171554.xlsx', 'Sheet', 'Motor X - Line 1','VariableNamingRule','preserve'));

% % Método 1 - FullTrap
% [espessura(1), Ue(1), H(1), delta_asterisco(1), theta(1), u{1}, ucl{1}, umais{1}, ymais{1}, uparedeclauser{1}] = calculateBLParametersFullTrap(pxraw1, pslocal(1), Cf(1), rho_air, nu_air);
% [espessura(2), Ue(2), H(2), delta_asterisco(2), theta(2), u{2}, ucl{2}, umais{2}, ymais{2}, uparedeclauser{2}] = calculateBLParametersFullTrap(pxraw2, pslocal(2), Cf(2), rho_air, nu_air);
% [espessura(3), Ue(3), H(3), delta_asterisco(3), theta(3), u{3}, ucl{3}, umais{3}, ymais{3}, uparedeclauser{3}] = calculateBLParametersFullTrap(pxraw3, pslocal(3), Cf(3), rho_air, nu_air);
% [espessura(4), Ue(4), H(4), delta_asterisco(4), theta(4), u{4}, ucl{4}, umais{4}, ymais{4}, uparedeclauser{4}] = calculateBLParametersFullTrap(pxraw4, pslocal(4), Cf(4), rho_air, nu_air);

% % Método 2 - TeoExp (Cf ábaco)
% [espessura(1), Ue(1), H(1), delta_asterisco(1), theta(1), u{1}, ucl{1}, umais{1}, ymais{1}, uparedeclauser{1}] = calculateBLParametersTeoExp(pxraw1, pslocal(1), Cfabaco(1), rho_air, nu_air);
% [espessura(2), Ue(2), H(2), delta_asterisco(2), theta(2), u{2}, ucl{2}, umais{2}, ymais{2}, uparedeclauser{2}] = calculateBLParametersTeoExp(pxraw2, pslocal(2), Cfabaco(2), rho_air, nu_air);
% [espessura(3), Ue(3), H(3), delta_asterisco(3), theta(3), u{3}, ucl{3}, umais{3}, ymais{3}, uparedeclauser{3}] = calculateBLParametersTeoExp(pxraw3, pslocal(3), Cfabaco(3), rho_air, nu_air);
% [espessura(4), Ue(4), H(4), delta_asterisco(4), theta(4), u{4}, ucl{4}, umais{4}, ymais{4}, uparedeclauser{4}] = calculateBLParametersTeoExp(pxraw4, pslocal(4), Cfabaco(4), rho_air, nu_air);

% Método 2 - TeoExp (Cf vonkarman1)
load("Cfvonkarman3.mat");
[espessura(1), Ue(1), H(1), delta_asterisco(1), theta(1), u{1}, ucl{1}, umais{1}, ymais{1}, uparedeclauser{1}, utau(1)] = calculateBLParametersTeoExp(pxraw1, pslocal(1), Cfvonkarman3(1), rho_air, nu_air);
[espessura(2), Ue(2), H(2), delta_asterisco(2), theta(2), u{2}, ucl{2}, umais{2}, ymais{2}, uparedeclauser{2}, utau(2)] = calculateBLParametersTeoExp(pxraw2, pslocal(2), Cfvonkarman3(2), rho_air, nu_air);
[espessura(3), Ue(3), H(3), delta_asterisco(3), theta(3), u{3}, ucl{3}, umais{3}, ymais{3}, uparedeclauser{3}, utau(3)] = calculateBLParametersTeoExp(pxraw3, pslocal(3), Cfvonkarman3(3), rho_air, nu_air);
[espessura(4), Ue(4), H(4), delta_asterisco(4), theta(4), u{4}, ucl{4}, umais{4}, ymais{4}, uparedeclauser{4}, utau(4)] = calculateBLParametersTeoExp(pxraw4, pslocal(4), Cfvonkarman3(4), rho_air, nu_air);


save('espessura', 'espessura');
save('Ue', 'Ue');
save('H', 'H');
save('delta_asterisco', 'delta_asterisco');
save('theta', 'theta');
save('u', 'u');
save('ucl', 'ucl');
save('umais', 'umais');
save('ymais', 'ymais');
save('uparedeclauser', 'uparedeclauser');
save('utau', 'utau');

%% vonkarman Cfs

load('x.mat');
load('theta.mat');
load('Ue.mat');
load('H.mat');

% % Finite differences
% diff_theta_x = diff(theta) ./ diff(x);
% dthetadx = [diff_theta_x(1), (diff_theta_x(1:end-1) + diff_theta_x(2:end))/2, diff_theta_x(end)];
% diff_Ue_x = diff(Ue) ./ diff(x);
% dUedx = [diff_Ue_x(1), (diff_Ue_x(1:end-1) + diff_Ue_x(2:end))/2, diff_Ue_x(end)];
% Cfvonkarman1 = 2.*(dthetadx + theta .* ((H + 2) ./ Ue) .* dUedx);
% save('Cfvonkarman1', 'Cfvonkarman1');

% % Splines
% theta_spline = spline(x, theta);
% theta_spline_derivative = fnder(theta_spline);
% dthetadx = ppval(theta_spline_derivative, x);
% Ue_spline = spline(x, Ue);
% Ue_spline_derivative = fnder(Ue_spline);
% dUedx = ppval(Ue_spline_derivative, x);
% Cfvonkarman2 = 2.*(dthetadx + theta .* ((H + 2) ./ Ue) .* dUedx);
% save('Cfvonkarman2', 'Cfvonkarman2');

% 2nd-degree polynomial fit
p_theta = polyfit(x, theta, 2);
p_Ue = polyfit(x, Ue, 2);
dthetadx_polyfit = 2 * p_theta(1) * x + p_theta(2);
dUedx_polyfit = 2 * p_Ue(1) * x + p_Ue(2);
Cfvonkarman3 = 2.*(dthetadx_polyfit + theta .* ((H + 2) ./ Ue) .* dUedx_polyfit);
save('Cfvonkarman3', 'Cfvonkarman3');


%% Plots configuration
resolution = 600;
markersize = 30;
%% I - Variação longitudinal do coeficiente de pressão estática

figure;
scatter(x./L1, Cp, markersize, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Scatter plot with blue filled circles
hold on; % Hold the plot to add more elements
p_adjust = polyfit(x./L1, Cp, 2); % Linear fit (adjustment line)
x_fit = linspace(min(x./L1), max(x./L1), 100); % Generate x values for fitted line
y_fit = polyval(p_adjust, x_fit); % Evaluate y values using the fitted line
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
xlabel('x/L', 'FontSize', 12); % Label for x-axis with font size
ylabel('C_p', 'FontSize', 12); % Label for y-axis with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 10); % Set font size for axis ticks
text(x./L1, Cp, compose(' %.3f', Cp), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
xlim([0, max(x./L1) + 0.1]); % Set x-axis limits based on data with padding
ylim([-0.4, 0.1]); % Set y-axis limits based on data with padding
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');
legend({'Dados experimentais'}, 'FontSize', 10, 'Location', 'best'); % Add legend

saveFigureAsPNG(resolution, 'I_Cp.png');

%% II - Variação longitudinal da velocidade exterior

figure;
scatter(x./L1, Uelocal, markersize, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Scatter plot with blue filled circles
hold on; % Hold the plot to add more elements
p_adjust = polyfit(x./L1, Uelocal, 2); % Linear fit (adjustment line)
x_fit = linspace(min(x./L1), max(x./L1), 100); % Generate x values for fitted line
y_fit = polyval(p_adjust, x_fit); % Evaluate y values using the fitted line
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
xlabel('x/L [-]', 'FontSize', 12); % Label for x-axis with font size
ylabel('U_e [m/s]', 'FontSize', 12); % Label for y-axis with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 10); % Set font size for axis ticks
text(x./L1, Uelocal, compose(' %.3f', Uelocal), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
xlim([0, max(x./L1) + 0.1]); % Set x-axis limits based on data with padding
ylim([20, 28]); % Set y-axis limits based on data with padding
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');
legend({'Dados experimentais'}, 'FontSize', 10, 'Location', 'best'); % Add legend

saveFigureAsPNG(resolution, 'II_Ue.png');

%% III - Região da camada da parede de um dos perfis de velocidade nas coordenadas de Clauser

perfil = 3; % alterar para ver outros
figure;
scatter(uparedeclauser{perfil}(:,1), uparedeclauser{perfil}(:,2), markersize, 'black', 'LineWidth', 1); % Scatter plot
hold on; % Hold the plot to add a line
p = polyfit(uparedeclauser{perfil}(:,1), uparedeclauser{perfil}(:,2), 2); % Second-degree polynomial fit
x_fit = linspace(min(uparedeclauser{perfil}(:,1)), max(uparedeclauser{perfil}(:,1)), 100); % Generate x values for fitting
y_fit = polyval(p, x_fit); % Evaluate the polynomial at x values
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
xlabel('Re_y', 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u/Ue', 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Helvetica', 'FontSize', 10); % Set font style and size for axis ticks
grid on; % Show gridlines
box on; % Show box around plot
ax = gca; % Get current axis
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w'); % Set figure background color to white
x_padding = 0.5 * (max(uparedeclauser{perfil}(:,1)) - min(uparedeclauser{perfil}(:,1)));
xlim([min(uparedeclauser{perfil}(:,1)) - x_padding, max(uparedeclauser{perfil}(:,1)) + x_padding]); % Set dynamic x-axis limits
y_padding = 0.8 * (max(uparedeclauser{perfil}(:,2)) - min(uparedeclauser{perfil}(:,2)));
ylim([min(uparedeclauser{perfil}(:,2)) - y_padding, max(uparedeclauser{perfil}(:,2)) + y_padding]); % Set dynamic y-axis limits
legend(['X', num2str(perfil), ' - Região da parede'], 'Location', 'best', 'FontSize', 10); % Add legend

saveFigureAsPNG(resolution, 'III_clauserparede.png');

%% IV - Os perfis de velocidade em escalas lineares

load('espessura.mat')
load('Ue.mat')
load('ucl.mat')

% Individuais

% X1
figure;
scatter(ucl{1}(:,1) / espessura(1), ucl{1}(:,2) / Ue(1), 15, 'filled', 'LineWidth', 0.5); % Customize scatter points
xlabel('y/\delta [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u/Ue [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');
xlim([0, 1]); % Set x-axis limits from 0 to 1
ylim([0.6, 1]); % Set y-axis limits from 0 to 1
legend('X1', 'Location', 'best', 'FontName', 'Arial','FontSize', 10);

saveFigureAsPNG(resolution, 'IV_X1_linear.png');

% X2
figure;
scatter(ucl{2}(:,1) / espessura(2), ucl{2}(:,2) / Ue(2), 15, 'filled', 'LineWidth', 0.5); % Customize scatter points
xlabel('y/\delta [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u/Ue [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');
xlim([0, 1]); % Set x-axis limits from 0 to 1
ylim([0.6, 1]); % Set y-axis limits from 0 to 1
legend('X2', 'Location', 'best', 'FontName', 'Arial','FontSize', 10);

saveFigureAsPNG(resolution, 'IV_X2_linear.png');

% X3
figure;
scatter(ucl{3}(:,1) / espessura(3), ucl{3}(:,2) / Ue(3), 15, 'filled', 'LineWidth', 0.5); % Customize scatter points
xlabel('y/\delta [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u/Ue [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');
xlim([0, 1]); % Set x-axis limits from 0 to 1
ylim([0.6, 1]); % Set y-axis limits from 0 to 1
legend('X3', 'Location', 'best', 'FontName', 'Arial','FontSize', 10);

saveFigureAsPNG(300, 'IV_X3_linear.png');

% X4
figure;
scatter(ucl{4}(:,1) / espessura(4), ucl{4}(:,2) / Ue(4), 15, 'filled', 'LineWidth', 0.5); % Customize scatter points
xlabel('y/\delta [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u/Ue [-]', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');
xlim([0, 1]); % Set x-axis limits from 0 to 1
ylim([0.6, 1]); % Set y-axis limits from 0 to 1
legend('X4', 'Location', 'best', 'FontName', 'Arial','FontSize', 10);

saveFigureAsPNG(resolution, 'IV_X4_linear.png');

% Conjunto
figure;
hold on; % Hold the plot for multiple data sets
for idx = 1:numel(ucl)
    scatter(ucl{idx}(:,1) / espessura(idx), ucl{idx}(:,2) / Ue(idx), 15, 'filled', 'LineWidth', 0.5);
end
xlabel('y/\delta', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u/Ue', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Arial', 'FontSize', 10); % Set font style and size for axis ticks
grid on; % Show gridlines
box on;
ax = gca;
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w'); % Set figure background color to white
xlim([0, 1 + 0.1]); % Set x-axis limits with padding
ylim([0.6, 1]); % Set y-axis limits with padding
legends = cellfun(@(x) sprintf('X%d', x), num2cell(1:numel(ucl)), 'UniformOutput', false);
legend(legends, 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);
hold off; % Release the plot hold

saveFigureAsPNG(resolution, 'IV_todos_lineares.png');

%% V - Variação longitudinal da espessura da camada limite, delta asterisco, teta e H

load('espessura.mat');
load('delta_asterisco.mat');
load('theta.mat');
load('H.mat');
load('x.mat');

% Espessura
xx = x./L1;
yy = espessura * 1000;
figure;
scatter(xx, yy, markersize,'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Scatter plot with blue filled circles
hold on; % Hold the plot to add more elements
p_adjust = polyfit(xx, yy, 2);
x_fit = linspace(min(xx), max(xx), 100); % Generate x values for fitted line
y_fit = polyval(p_adjust, x_fit); % Evaluate y values using the fitted line
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
xlabel('x/L', 'FontSize', 12); % Label for x-axis with font size
ylabel('\delta [mm]', 'FontSize', 12); % Label for y-axis with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 10); % Set font size for axis ticks
text(x./L1, yy, compose(' %.3f', yy), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
xlim([0, 1]); % Set x-axis limits based on data with padding
ylim([0, 20]); % Set y-axis limits based on data with padding
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');

saveFigureAsPNG(resolution, 'V_espessura.png');

% Espessura de deslocamento
xx = x./L1;
yy = delta_asterisco * 1000;
figure;
scatter(xx, yy, markersize,'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Scatter plot with blue filled circles
hold on; % Hold the plot to add more elements
p_adjust = polyfit(xx, yy, 2);
x_fit = linspace(min(xx), max(xx), 100); % Generate x values for fitted line
y_fit = polyval(p_adjust, x_fit); % Evaluate y values using the fitted line
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
xlabel('x/L', 'FontSize', 12); % Label for x-axis with font size
ylabel('\delta* [mm]', 'FontSize', 12); % Label for y-axis with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 10); % Set font size for axis ticks
text(x./L1, yy, compose(' %.3f', yy), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
xlim([0, 1]); % Set x-axis limits based on data with padding
ylim([0, 5]); % Set y-axis limits based on data with padding
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');

saveFigureAsPNG(resolution, 'V_delta_asterisco.png');

% Espessura de défice de quantidade de movimento
xx = x./L1;
yy = theta * 1000;
figure;
scatter(xx, yy, markersize,'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Scatter plot with blue filled circles
hold on; % Hold the plot to add more elements
p_adjust = polyfit(xx, yy, 2);
x_fit = linspace(min(xx), max(xx), 100); % Generate x values for fitted line
y_fit = polyval(p_adjust, x_fit); % Evaluate y values using the fitted line
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
xlabel('x/L ', 'FontSize', 12); % Label for x-axis with font size
ylabel('\theta [mm]', 'FontSize', 12); % Label for y-axis with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 10); % Set font size for axis ticks
text(x./L1, yy, compose(' %.3f', yy), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
xlim([0, 1]); % Set x-axis limits based on data with padding
ylim([0, 3]); % Set y-axis limits based on data with padding
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');

saveFigureAsPNG(resolution, 'V_theta.png');

%% VI - Variação longitudinal do fator de forma da camada limite 
% Fator de forma
xx = x./L1;
yy = H;
figure;
scatter(xx, yy, markersize,'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b'); % Scatter plot with blue filled circles
hold on; % Hold the plot to add more elements
p_adjust = polyfit(xx, yy, 2);
x_fit = linspace(min(xx), max(xx), 100); % Generate x values for fitted line
y_fit = polyval(p_adjust, x_fit); % Evaluate y values using the fitted line
plot(x_fit, y_fit, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]); % Dotted line with transparency
xlabel('x/L', 'FontSize', 12); % Label for x-axis with font size
ylabel('H', 'FontSize', 12); % Label for y-axis with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 10); % Set font size for axis ticks
text(x./L1, yy, compose(' %.3f', yy), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
xlim([0, 1]); % Set x-axis limits based on data with padding
ylim([0, 3]); % Set y-axis limits based on data with padding
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w');

saveFigureAsPNG(resolution, 'VI_H.png');

%% VII - Variação longitudinal dos valores do coeficiente de tensão de corte superficial

load('x.mat');
load('Cfabaco.mat');
load('Cfvonkarman3.mat');

xx = x ./ L1;
yy1 = Cfabaco * 1000;
yy2 = Cfvonkarman3 * 1000;
figure;
scatter(xx, yy1, markersize, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
hold on;
p_adjust1 = polyfit(xx, yy1, 2);
x_fit1 = linspace(min(xx), max(xx), 100);
y_fit1 = polyval(p_adjust1, x_fit1);
plot(x_fit1, y_fit1, 'k--', 'LineWidth', 0.5, 'Color', [0 0 0 0.5]);
scatter(xx, yy2, markersize, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
p_adjust2 = polyfit(xx, yy2, 2);
x_fit2 = linspace(min(xx), max(xx), 100);
y_fit2 = polyval(p_adjust2, x_fit2);
plot(x_fit2, y_fit2, 'r--', 'LineWidth', 0.5, 'Color', [1 0 0 0.5]);
xlabel('x/L', 'FontSize', 12);
ylabel('C_f \times 10^3', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 10);
text(xx, yy1, compose(' %.2f', yy1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'b');
text(xx, yy2, compose(' %.2f', yy2), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'r');
xlim([0, 1]);
ylim([0, 5]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.GridColor = 'k';
ax.GridAlpha = 0.3;
ax.GridLineStyle = ':';
set(gcf, 'Color', 'w');
legend({'Ábaco de Clauser', '', 'von Kármán', ''}, 'Location', 'Best', 'FontSize', 10);
hold off;

saveFigureAsPNG(resolution, 'VII_Cf.png');

% por balanço dos diversos termos que figuram na equação integral de von-Kárman
load('x.mat');
load('Cfabaco.mat');
load('Cfvonkarman1.mat');
load('Cfvonkarman2.mat');
load('Cfvonkarman3.mat');

figure;
hold on;
plot(x, Cfabaco, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cfabaco');
plot(x, Cfvonkarman1, 's--', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cfvonkarman1');
plot(x, Cfvonkarman2, 'd-.', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cfvonkarman2');
plot(x, Cfvonkarman3, 's-.', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Cfvonkarman3');
xlabel('Index', 'FontSize', 12); % Label for x-axis with font size
ylabel('Value', 'FontSize', 12); % Label for y-axis with font size
legend('Location', 'Best', 'FontSize', 10); % Legend with font size
grid on; % Show grid
box on; % Show box around plot
set(gca, 'FontSize', 12); % Set font size for axis ticks
ax = gca; % Get current axis
ax.XColor = 'k'; % Set x-axis color to black
ax.YColor = 'k'; % Set y-axis color to black
ax.GridColor = [0.1, 0.1, 0.1]; % Set grid color
ax.GridAlpha = 0.5; % Set grid transparency
ax.GridLineStyle = '--'; % Set grid line style
ax.XGrid = 'on'; % Turn on x-axis grid
ax.YGrid = 'on'; % Turn on y-axis grid
set(gcf, 'Color', 'w');
hold off;

saveFigureAsPNG(resolution, 'Cfs_comparison.png');

%% VIII - Perfis semi-logarítmicos

load('umais.mat');
load('ymais.mat');
load('utau');

% Individual X1

kappa = 0.41; % von Karman constant
B = 5.2; % constant for the law of the wall
law_of_wall = @(y_plus) 1/kappa * log(y_plus) + B;

% X1
figure;
hold on;
scatter(ymais{1}(:), umais{1}(:), 15, 'filled', 'LineWidth', 0.5);
y_plus_range = logspace(0.1, 3); % Define range for log scale
u_plus = law_of_wall(y_plus_range);
plot(y_plus_range, u_plus, 'k--', 'LineWidth', 1.5);
xlabel('y^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('u^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k';
ax.GridAlpha = 0.3;
ax.GridLineStyle = ':';
xlim([30, 1000]);
%ylim([0.6, 1]);
set(gcf, 'Color', 'w');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend('Perfil X1', 'Lei da parede', 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);
hold off;

saveFigureAsPNG(resolution, 'VIII_X1semilog.png');

% X2
figure;
hold on;
scatter(ymais{2}(:), umais{2}(:), 15, 'filled', 'LineWidth', 0.5);
y_plus_range = logspace(0.1, 3); % Define range for log scale
u_plus = law_of_wall(y_plus_range);
plot(y_plus_range, u_plus, 'k--', 'LineWidth', 1.5);
xlabel('y^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('u^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k';
ax.GridAlpha = 0.3;
ax.GridLineStyle = ':';
xlim([30, 1000]);
%ylim([0.6, 1]);
set(gcf, 'Color', 'w');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend('Perfil X1', 'Lei da parede', 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);
hold off;

saveFigureAsPNG(resolution, 'VIII_X2semilog.png');

% X3
figure;
hold on;
scatter(ymais{3}(:), umais{3}(:), 15, 'filled', 'LineWidth', 0.5);
y_plus_range = logspace(0.1, 3); % Define range for log scale
u_plus = law_of_wall(y_plus_range);
plot(y_plus_range, u_plus, 'k--', 'LineWidth', 1.5);
xlabel('y^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('u^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k';
ax.GridAlpha = 0.3;
ax.GridLineStyle = ':';
xlim([30, 1000]);
%ylim([0.6, 1]);
set(gcf, 'Color', 'w');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend('Perfil X1', 'Lei da parede', 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);
hold off;

saveFigureAsPNG(resolution, 'VIII_X3semilog.png');

% X4
figure;
hold on;
scatter(ymais{4}(:), umais{4}(:), 15, 'filled', 'LineWidth', 0.5);
y_plus_range = logspace(0.1, 3); % Define range for log scale
u_plus = law_of_wall(y_plus_range);
plot(y_plus_range, u_plus, 'k--', 'LineWidth', 1.5);
xlabel('y^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('u^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
grid on;
box on;
ax = gca;
ax.GridColor = 'k';
ax.GridAlpha = 0.3;
ax.GridLineStyle = ':';
xlim([30, 1000]);
%ylim([0.6, 1]);
set(gcf, 'Color', 'w');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend('Perfil X1', 'Lei da parede', 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);
hold off;

saveFigureAsPNG(resolution, 'VIII_X4semilog.png');


% Conjunto
figure;
hold on; % Hold the plot for multiple data sets
for idx = 1:numel(ymais)
    scatter(ymais{idx}(:), umais{idx}(:), 15, 'filled', 'LineWidth', 0.5);
end
y_plus_range = logspace(0.1, 3); % Define range for log scale
u_plus = law_of_wall(y_plus_range);
plot(y_plus_range, u_plus, 'k--', 'LineWidth', 1.5);
xlabel('y^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for x-axis
ylabel('u^+', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal'); % Label for y-axis
set(gca, 'FontName', 'Arial', 'FontSize', 10); % Set font style and size for axis ticks
grid on; % Show gridlines
box on;
ax = gca;
ax.GridColor = 'k'; % Set grid color to black
ax.GridAlpha = 0.3; % Set grid transparency
ax.GridLineStyle = ':'; % Set grid line style
set(gcf, 'Color', 'w'); % Set figure background color to white
xlim([40, 1100]);
% ylim([0.6, 1]);
set(gca, 'XScale', 'log');
legends = cellfun(@(x) sprintf('X%d', x), num2cell(1:numel(ucl)), 'UniformOutput', false);
legend([legends, 'Lei da parede'], 'Location', 'best', 'FontName', 'Arial', 'FontSize', 10);
hold off; % Release the plot hold

saveFigureAsPNG(resolution, 'VIII_todossemilog.png');

%% outros

load('delta_asterisco.mat');
load('theta.mat');
load('H.mat');
load('Cfabaco.mat');
load('Cfvonkarman1.mat');
load('Cfvonkarman2.mat');
load('Cfvonkarman3.mat');

paramnames = {'delta_asterisco', 'theta', 'H', 'Cf_abaco', 'Cf_vonkarman1', 'Cf_vonkarman2', 'Cf_vonkarman3'};
stations = {'X1', 'X2', 'X3', 'X4'};
param = [delta_asterisco; theta; H; Cfabaco; Cfvonkarman1; Cfvonkarman2; Cfvonkarman3];
paramintTable = array2table(param, 'RowNames', paramnames, 'VariableNames', stations);
disp(paramintTable);