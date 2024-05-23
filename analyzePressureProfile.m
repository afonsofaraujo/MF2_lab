function analyzePressureProfile(p_X, p0, rho_air, nu_air)
    % retirar pontos iniciais
    pmax = max(p_X(:,2));
    pmin = min(p_X(:,2));
    deltamax = pmax - pmin;

    i = 1;
    while i <= length(p_X(:,2)) && (p_X(i,2) - pmin) <= 0.01 * deltamax
        i = i + 1;
    end
    
    if i > 1
        p_X(:,1) = p_X(:,1) - p_X(i-1,1);
        p_X = p_X(i-1:end,:);
        fprintf('Foram retirados os %i primeiros pontos\n', i-1);
    else
        % If i is 1 and the condition is not met for any point, handle it here
        % For example, you might choose to not remove any points.
        fprintf('Nenhum ponto removido\n');
    end

    % correção das leituras dos deslocamentos verticais
    McMillan = (0.5 + 0.15) * 1;
    p_X(:,1) = p_X(:,1) + McMillan;

    % passar a velocidades e y a [m]
    u_X = [p_X(:,1) / 1000, sqrt(2*(p_X(:,2)-p0) / rho_air)];

    % definir velocidade exterior como velocidade máxima e fim da camada limite
    Ue = max(u_X(:,2));
    
    % Find the index where the velocity drops below 0.985 * Ue
    i = length(u_X(:,2));
    while i >= 1 && u_X(i,2) >= 0.985 * Ue
        i = i - 1;
    end
    
    % Remove points beyond the detected index
    if i < length(u_X(:,2))
        fprintf('Foram retirados os %i últimos pontos\n', length(u_X(:,2)) - i);
        u_cl_X = u_X(1:i,:);
    else
        fprintf('Nenhum ponto removido\n');
        u_cl_X = u_X; % Optionally keep the original u_X if no removal occurs
    end

    % coordenadas de Clauser
    u_clauser_X = [Ue * u_cl_X(:,1) / nu_air, u_cl_X(:,2) / Ue];

    % filtrar cordenadas de clauser <15% da CL
    espessura = u_cl_X(end,1);
    u_cl_X_filtered = u_cl_X(u_cl_X(:,1) < 0.15 * espessura, :);
    u_clauser_X_filtered = [Ue * u_cl_X_filtered(:,1) / nu_air, u_cl_X_filtered(:,2) / Ue];

    % Cf com a interpolação do ábaco de Clauser
    Cf = 0.004;

    % adimensionalizar
    u_tau = Ue * sqrt(Cf/2);
    umais = u_cl_X(:,2) / u_tau;
    ymais = (u_cl_X(:,1) * u_tau) / nu_air;
    deltamais = (espessura * u_tau) / nu_air;

    % adicionar o ponto (0,0)
    umais = [0; umais];
    ymais = [0; ymais];

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

    % Display or return relevant results
    disp(H);
end
