function [espessura, Ue, H, delta_asterisco, teta, uX, uXcl, umais, ymais] = calculateBLParametersFullTrap(pXraw, ps, Cf, rho_air, nu_air)
    % clean pX
    pX = [pXraw(:,1) / 1000, pXraw(:,2)]; % passar [m]
    i = 1;
    while i <= length(pX(:,2)) && (pX(i,2) - min(pX(:,2))) <= 0.015 * (max(pX(:,2)) - min(pX(:,2)))
        i = i + 1;
    end
    if i > 1
        pX(:,1) = pX(:,1) - pX(i-1,1);
        pX = pX(i-1:end,:);
        fprintf('Foram retirados os %i primeiros pontos\n', i-1);
    else
        fprintf('Nenhum ponto removido\n');
    end
    % correção das leituras dos deslocamentos verticais
    McMillan = (0.0005 + 0.00015) * 1;
    pX(:,1) = pX(:,1) + McMillan;
    % velocidades
    uX = [ pX(:,1), sqrt( 2 * (pX(:,2) - ps) / rho_air)];
    Ue = max(uX(:,2));
    i = length(uX(:,2));
    while i >= 1 && uX(i,2) >= 0.985 * Ue
        i = i - 1;
    end
    if i < length(uX(:,2))
        fprintf('Foram retirados os %i últimos pontos\n', length(uX(:,2)) - i);
        uXcl = uX(1:i,:);
    else
        fprintf('Nenhum ponto removido\n');
        uXcl = uX; % Optionally keep the original u_X if no removal occurs
    end
    espessura = uXcl(end,1);

    % Adimensionalizar
    u_tau = Ue * sqrt(Cf / 2);
    umais = uXcl(:,2) / u_tau;
    ymais = (uXcl(:,1) * u_tau) / nu_air;
    deltamais = (espessura * u_tau) / nu_air;

    % Adicionar o ponto (0,0)
    umais = [0; umais];
    ymais = [0; ymais];

    % Espessura de deslocamento - camada interior/exterior
    delta_asterisco = espessura - (nu_air / Ue) * int_trap(ymais, umais, 0, deltamais, 0);

    % Espessura de défice de quantidade de movimento - camada interior/exterior
    teta = (nu_air / Ue) * int_trap(ymais, umais, 0, deltamais, 0) - (nu_air / Ue) * sqrt(Cf / 2) * int_trap(ymais, umais.^2, 0, deltamais, 0);

    % Fator de forma
    H = delta_asterisco / teta;
end