clear
clc

% Funções:
function cdf = finitecdf(s, p, x)
    % Variável aleatória finita X:
    % vetor s com o espaço amostral
    % elementos {s(1), s(2), ...}
    % vetor p com as probabilidades
    % p(i) = P[X = s(i)]
    % A saída é o vetor cdf:
    % cdf(i) = P[X <= x(i)]

    cdf = zeros(size(x));
    for i = 1:length(x)
        pxi = sum(p(s <= x(i)));
        cdf(i) = pxi;
    end
end

function rho = finitecoeff(SX, SY, PXY)
    % Uso: rho = finitecoeff(SX, SY, PXY)
    % Calcula o coeficiente de correlação rho de
    % variáveis aleatórias finitas X e Y
    ex = finiteexp(SX, PXY);
    vx = finitevar(SX, PXY);
    ey = finiteexp(SY, PXY);
    vy = finitevar(SY, PXY);
    R = finiteexp(SX .* SY, PXY);
    rho = (R - ex * ey) / sqrt(vx * vy);
end

function covxy = finitecov(SX, SY, PXY)
    % Uso: covxy = finitecov(SX, SY, PXY)
    % Retorna a covariância das variáveis aleatórias finitas X e Y
    % dadas pelas grades SX, SY e PXY
    
    ex = finiteexp(SX, PXY);
    ey = finiteexp(SY, PXY);
    R = finiteexp(SX .* SY, PXY);
    covxy = R - ex * ey;
end

function ex = finiteexp(sx, px)
    % Uso: ex = finiteexp(sx, px)
    % Retorna o valor esperado E[X]
    % da variável aleatória finita X descrita
    % pelos valores amostrais sx e probabilidades px
    ex = sum((sx(:)) .* (px(:)));
end
  
function pmf = finitepmf(sx, px, x)
    % Variável aleatória finita X:
    % vetor sx com o espaço amostral
    % elementos {sx(1), sx(2), ...}
    % vetor px com as probabilidades
    % px(i) = P[X = sx(i)]
    % A saída é o vetor pmf:
    % pmf(i) = P[X = x(i)]

    pmf = zeros(size(x(:)));
    for i = 1:length(x)
        pmf(i) = sum(px(sx == x(i)));
    end
end

function x = finiterv(s, p, m)
    % Retorna m amostras
    % de uma variável aleatória finita (s, p)
    % s = s(:); p = p(:);
    r = rand(m, 1);
    x = zeros(m, 1);

    for i = 1:m
        for j = 1:length(s)
            if r(i) <= sum(p(1:j))
                x(i) = s(j);
                break;
            end
        end
    end
end
  
function v = finitevar(sx, px)
    % Uso: v = finitevar(sx, px)
    % Retorna a variância Var[X]
    % de variáveis aleatórias finitas X descritas por
    % amostras sx e probabilidades px
    ex2 = finiteexp(sx.^2, px);
    ex = finiteexp(sx, px);
    v = ex2 - (ex^2);
end

%Input:
% Definindo a grade de amostras e probabilidades
sx = [1, 2, 3, 4]; % Exemplo de espaço amostral
px = [0.1, 0.2, 0.3, 0.4]; % Exemplo de probabilidades
m = input('Digite o numero de amostras: '); % Número de amostras

% Gerando amostras da variável aleatória finita
x_values = finiterv(sx, px, m);


% Output:
% Plotando a PMF
subplot(2, 1, 1);
pmf = finitepmf(sx, px, sx); % Usei sx como valores de x para a PMF
stem(sx, pmf, 'Marker', 'o', 'LineStyle', 'none');
title('Função Massa de Probabilidade (PMF)');
xlabel('X');
ylabel('Probabilidade');
grid on;

% Plotando a CDF
subplot(2, 1, 2);
cdf = finitecdf(sx, px, sx); % Usei sx como valores de x para a CDF
stairs(sx, cdf, 'Marker', 'o', 'LineStyle', '-');
title('Função de Distribuição Cumulativa (CDF)');
xlabel('X');
ylabel('Probabilidade acumulada');
grid on;
