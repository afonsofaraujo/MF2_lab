function [espessura, Ue, H, delta_asterisco, theta, u, ucl, umais, ymais, uparedeclauser, utau] = calculateBLParametersTeoExp(praw, ps, Cf, rho, nu)
    

    % clean pX
    p = [praw(:,1) / 1000, praw(:,2)]; % passar [m]
    i = 1;
    while i <= length(p(:,2)) && (p(i,2) - min(p(:,2))) <= 0.015 * (max(p(:,2)) - min(p(:,2)))
        i = i + 1;
    end
    if i > 1
        p(:,1) = p(:,1) - p(i-1,1);
        p = p(i-1:end,:);
        fprintf('Foram retirados os %i primeiros pontos\n', i-1);
    else
        fprintf('Nenhum ponto removido\n');
    end

    % correção das leituras dos deslocamentos verticais
    McMillan = (0.0005 + 0.00015) * 1;
    p(:,1) = p(:,1) + McMillan;
    
    % velocidades
    u = [ p(:,1), sqrt( 2 * (p(:,2) - ps) / rho)];
    Ue = max(u(:,2));

    i = length(u(:,2));
    while i >= 1 && u(i,2) >= 0.985 * Ue
        i = i - 1;
    end
    if i < length(u(:,2))
        fprintf('Foram retirados os %i últimos pontos\n', length(u(:,2)) - i);
        ucl = u(1:i,:);
    else
        fprintf('Nenhum ponto removido\n');
        ucl = u; % Optionally keep the original u_X if no removal occurs
    end
    espessura = ucl(end,1);

    % Adimensionalizar
    utau = Ue * sqrt(Cf / 2);
    umais = ucl(:,2) / utau;
    ymais = (ucl(:,1) * utau) / nu;
    deltamais = (espessura * utau) / nu;

    % filtrar cordenadas de clauser y < 15% da CL e y+ > 50
    uparede = ucl(ucl(:,1) < 0.15 * espessura & ucl(:,1) > 50 * nu / utau, :);
    uparedeclauser = [Ue*uparede(:,1) / nu, uparede(:,2) / Ue];

    % constantes de Von-Karman
    k = 0.41;
    C = 0.52;

    %% delta_asterisco
    % coles y = [0, (50 * nu) / u_tau]
    intcoles1 = 540.6;
    coles1 = (50 * nu) / utau - (nu / Ue) * intcoles1;
    % parede y = [(50 * nu) / u_tau, 0.15 * deltamais * nu / u_tau]
    intparede1 = (1 / k) * (0.15 * deltamais * log(0.15 * deltamais) - 50 * log(50) - 0.15 * deltamais + 50) + C * (0.15 * deltamais - 50);
    parede1 = (0.15 * espessura - (50 * nu) / utau) - (nu / Ue) * intparede1;
    % exterior y = [0.15*espessura, espessura]
    exterior1 = int_trap(ucl(:,1), 1 - ucl(:,2) ./ Ue, 0.15*espessura, espessura, 0);

    delta_asterisco = coles1 + parede1 + exterior1;
    %% theta
    % coles y = [0, (50 * nu) / u_tau]
    intcoles2 = 6546;
    coles2 = (nu / Ue) * (intcoles1 - intcoles2 * sqrt(Cf / 2));
    % parede y = [(50 * nu) / u_tau, 0.15 * deltamais * nu / u_tau]
    a = (50 * nu) / utau;
    b = 0.15 * deltamais;
    intparede2 = (1 / k^2) * (2*(b - a) + C*k*(C * k - 2)*(b - a) + b*log(b)*(log(b) + 2*C*k - 2) - a*log(a)*(log(a) + 2*C*k - 2));
    parede2 = (nu / Ue) * intparede1 - (nu / Ue) * sqrt(Cf / 2) * intparede2;
    % exterior y = [0.15 * espessura, espessura]
    exterior2 = int_trap(ucl(:,1), (ucl(:,2) / Ue) - (ucl(:,2) / Ue).^2, 0.15*espessura, espessura, 0);

    theta = coles2 + parede2 + exterior2;
    %% Fator de forma
    H = delta_asterisco / theta;
end