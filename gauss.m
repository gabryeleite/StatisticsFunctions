clear
clc

% Funções:
function f = gausspdf(mu, sigma, x)
    f = exp(-(x - mu).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);
end

function f = gausscdf(mu, sigma, x)
  f = 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));
end

function x = gaussrv(mu, sigma, m)
    x = mu + sigma * randn(m, 1);
end

function x = gaussvector(mu, C, m)
    % Saída: m vetores gauss, cada um com média mu e matriz de covariância C

    if (min(size(C)) == 1)
        C = toeplitz(C);
    end

    n = size(C, 2);

    if (length(mu) == 1)
        mu = mu * ones(n, 1);
    end

    [U, D, V] = svd(C);
    x = V * (D^(0.5)) * randn(n, m) + (mu(:) * ones(1, m));
end

function f = gaussvectorpdf(mu, C, a)
    n = length(a);
    z = a(:) - mu(:);
    f = exp(-z' * inv(C) * z) / sqrt((2 * pi)^n * det(C));
end

% Input:
mu = input('Informe o parâmetro mu (inteiro positivo e dentro dos dados fornecidos): ');
C = sigma = input('Informe a probabilidade sigma (entre 0 e 1): ');
%a = m = input('Informe o número de amostras (m): ');
% Quando sigma é muito pequeno, a função de densidade de probabilidade (PDF) se torna muito "pico",
% e quando sigma é grande, a PDF se espalha. (Podemos conferir nos gráficos)
fid = fopen('dados.txt','r');
x = fscanf(fid, '%d');
fclose(fid);
%x = gaussrv(mu, sigma, m)
%y = gaussvector(mu, C, m)
valores_pdf = gausspdf(mu, sigma, x);
valores_cdf = gausscdf(mu, sigma, x);
%valores_vet = gaussvectorpdf(mu, C, a)

% Output:
% GaussPDF
subplot(2,1,1);
stem(x, valores_pdf, 'LineWidth', 2);
title('Gauss PDF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
xticks(x);
grid on;

% GausslCDF
subplot(2,1,2);
stairs(x, valores_cdf, 'LineWidth', 2); % stairs melhor representação para CDF, como se trata de uma soma acumulada
title('Gauss CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
xticks(x);
grid on;
