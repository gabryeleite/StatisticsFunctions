clear
clc

function f = uniformpdf(a, b, x)
  % Uso: f = uniformpdf(a, b, x)
  % Retorna a PDF de uma variável aleatória uniforme contínua avaliada em x
  f = ((x >= a) & (x < b)) / (b - a);
end

function F = uniformcdf(a, b, x)
  % Uso: F = uniformcdf(a, b, x)
  % Retorna a CDF de uma variável aleatória uniforme contínua avaliada em x
  F = x .* ((x >= a) & (x < b)) / (b - a);
  F = F + 1.0 * (x >= b);
  F(F > 1) = 1;
end

function x = uniformrv(a, b, m)
  % Uso: x = uniformrv(a, b, m)
  % Retorna m amostras de uma variável aleatória uniforme (a, b)
  x = a + (b - a) * rand(m, 1);
end

% Input:
a = input('Digite o limite inferior: '); % Limite inferior da distribuição
b = input('Digite o limite superior: '); % Limite superior da distribuição
%m = input('Digite o numero de amostras: ');
%x = uniformrv(a, b, m);
x = a : 0.1 : b;

%fid = fopen('dados.txt','r');
%x = fscanf(fid, '%d');
%fclose(fid);

% Calcular PDF e CDF
valores_pdf = uniformpdf(a, b, x);
valores_cdf = uniformcdf(a, b, x);

% Output:
% UniformePDF
subplot(2, 1, 1);
plot(x, valores_pdf, 'LineWidth', 2);
title('Uniforme PDF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
%xticks(x);
grid on;

% UniformeCDF
subplot(2, 1, 2);
plot(x, valores_cdf, 'LineWidth', 2);
title('Uniforme CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
%xticks(x);
grid on;
