close all
clear

% propriedades do ar (21ºC)
pT = 100658;
rho_air = 1.201;
nu_air = 1.515e-5;
mu_air = 1.819e-5;
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
x = [0.225; 0.475; 0.725; 0.975];
l = latm - manometro;
% ponto de referência dentro do ventilador
pTref = (latm - lpTref) * rho_fm * g * sind(beta);
psref = (latm - lpsref) * rho_fm * g * sind(beta);
qref = pTref - psref;
Ueref = sqrt(2 * qref / rho_air);
pslocal = l * rho_fm * g * sind(beta);
qlocal = pTref - pslocal;
Uelocal = sqrt(2 * qlocal / rho_air);
Cp = abs((pslocal - psref) / qref);





%% clean px
px = [pxraw(:,1) / 1000, pxraw(:,2)]; % passar [m]
i = 1;
while i <= length(px(:,2)) && (px(i,2) - min(px(:,2))) <= 0.015 * (max(px(:,2)) - min(px(:,2)))
    i = i + 1;
end
if i > 1
    px(:,1) = px(:,1) - px(i-1,1);
    px = px(i-1:end,:);
    fprintf('Foram retirados os %i primeiros pontos\n', i-1);
else
    fprintf('Nenhum ponto removido\n');
end
%% correção das leituras dos deslocamentos verticais
McMillan = (0.0005 + 0.00015) * 1;
px(:,1) = px(:,1) + McMillan;
%% velocidades
ux = [ px(:,1), sqrt( 2 * (px(:,2) - pslocal(1)) / rho_air)];
Ue = max(ux(:,2));
i = length(ux(:,2));
while i >= 1 && ux(i,2) >= 0.985 * Ue
    i = i - 1;
end
if i < length(ux(:,2))
    fprintf('Foram retirados os %i últimos pontos\n', length(ux(:,2)) - i);
    uxcl = ux(1:i,:);
else
    fprintf('Nenhum ponto removido\n');
    uxcl = ux; % Optionally keep the original u_x if no removal occurs
end
espessura = uxcl(end,1);

% coordenadas de Clauser
uxclauser = [Ue * uxcl(:,1) / nu_air, uxcl(:,2) / Ue];

% Cf com a interpolação do ábaco de Clauser
Cf = 0.0039;

% adimensionalizar
u_tau = Ue * sqrt(Cf/2);
umais = uxcl(:,2) / u_tau;
ymais = (uxcl(:,1) * u_tau) / nu_air;
deltamais = (espessura * u_tau) / nu_air;

% adicionar o ponto (0,0)
umais = [0; umais];
ymais = [0; ymais];

% filtrar cordenadas de clauser y < 15% da CL e y+ > 50
uxparede = uxcl(uxcl(:,1) < 0.15 * espessura & uxcl(:,1) > 50 * nu_air / u_tau, :);
uxparedeclauser = [Ue*uxparede(:,1) / nu_air, uxparede(:,2) / Ue];

% Format the values into a comma-separated string
Re_y_string = sprintf('%.15g, ', uparedeclauser{3}(:,1));
Re_y_string = Re_y_string(1:end-2);
u_Ue_string = sprintf('%.15g, ', uparedeclauser{3}(:,2));
u_Ue_string = u_Ue_string(1:end-2);
disp(Re_y_string);
disp(u_Ue_string);

% espessura de deslocamento - camada interior/exterior
delta_asterisco_1 = 0.15 * espessura - (nu_air / Ue) * int_trap(ymais, umais, 0, 0.15 * deltamais, 0);
delta_asterisco_2 = 0.85 * espessura - (nu_air / Ue) * int_trap(ymais, umais, 0.15 * deltamais, deltamais, 0);
delta_asterisco = espessura - (nu_air / Ue) * int_trap(ymais, umais, 0, deltamais, 0);

% subcamadas da camada interior adimensionais - teórico
delta_asterisco_coles = int_trap(ymais, umais, 0, 50, 0);
delta_asterisco_parede = int_trap(ymais, umais, 50, 0.15 * deltamais, 0);

% espessura de défice de quantidade de movimento - camada interior/exterior
teta_1 = (nu_air / Ue) * int_trap(ymais, umais, 0, 0.15 * deltamais, 0) - (nu_air / Ue) * sqrt(Cf / 2) * int_trap(ymais, umais.^2, 0, 0.15 * deltamais, 0);
teta_2 = (nu_air / Ue) * int_trap(ymais, umais, 0.15 * deltamais, deltamais, 0) - (nu_air / Ue) * sqrt(Cf / 2) * int_trap(ymais, umais.^2,  0.15 * deltamais, deltamais, 0);
teta = (nu_air / Ue) * int_trap(ymais, umais, 0, deltamais, 0) - (nu_air / Ue) * sqrt(Cf / 2) * int_trap(ymais, umais.^2, 0, deltamais, 0);

% fator de forma
H = delta_asterisco / teta;

disp(delta_asterisco);
disp(teta);
disp(H);

% plots

figure;

subplot(2, 3, 1);
scatter(px(:,1),px(:,2),"red");
ylabel("p [Pa]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
xlabel("y [mm]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
grid on

subplot(2, 3, 2);
scatter(ux(:,1),ux(:,2));
ylabel("v [m/s]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
xlabel("y[mm]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
grid on

subplot (2, 3, 3);
scatter(uxcl(:,1),uxcl(:,2));
ylabel("v [m/s]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
xlabel("y[mm]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
grid on

subplot(2, 3, 4);
scatter( uxparedeclauser(:,1), uxparedeclauser(:,2),"black","LineWidth",1);
set(gca, 'XScale', 'log');
xlabel("Re_y [-]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
ylabel("u/Ue [-]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal");
grid on

subplot(2, 3, 5);
scatter( uxclauser(:,1), uxclauser(:,2),"black","LineWidth",1)
set(gca, 'XScale', 'log');
xlabel("Re_y [-]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal")
ylabel("u/Ue [-]", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal")
grid on
