clear
clc

function f = gausspdf(mu, sigma, x)
  f = exp(-(x - mu).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);
end

function f = gausscdf(mu, sigma, x)
  f = 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));
end

function x = gaussrv(mu, sigma, m)
  x = mu + (sigma * randn(m, 1));
end

% Input:
a = input('Digite o limite inferior: '); % Limite inferior da distribuição
b = input('Digite o limite superior: '); % Limite superior da distribuição
mu = input('Informe o parâmetro mu (inteiro positivo): ');
sigma = input('Informe a probabilidade sigma (entre 0 e 1): ');
%m = input('Informe o número de amostras (m): ');
%x = gaussrv(mu, sigma, m);
x = a : 0.1 : b;

%fid = fopen('dados.txt','r');
%x = fscanf(fid, '%d');
%fclose(fid);

% Calcular PDF e CDF
valores_pdf = gausspdf(mu, sigma, x);
valores_cdf = gausscdf(mu, sigma, x);

% Output:
% GaussPDF
subplot(2,1,1);
plot(x, valores_pdf, 'LineWidth', 2);
title('Gauss PDF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
%xticks(x);
grid on;

% GausslCDF
subplot(2,1,2);
plot(x, valores_cdf, 'LineWidth', 2);
title('Gauss CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
%xticks(x);
grid on;
