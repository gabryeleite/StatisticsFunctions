clear
clc

function f = exponencialpdf(lambda, x)
  f = lambda * exp(-lambda * x);
  f = f .* (x >= 0);
end

function F = exponencialcdf(lambda, x)
  F = 1.0 - exp(-lambda * x);
end

function x = exponencialrv(lambda, m)
  x = -(1/lambda) * log(1 - rand(m, 1));
end

% Input:
a = input('Digite o limite inferior: '); % Limite inferior da distribuição
b = input('Digite o limite superior: '); % Limite superior da distribuição
lambda = input('Digite a taxa de chegada lambda (entre 0 e 1): ');
%m = input('Digite o numero de amostras: ');
%x = exponencialrv(lambda, m);
x = a : 0.1 : b;

%fid = fopen('dados.txt','r');
%x = fscanf(fid, '%d');
%fclose(fid);

Epdf = exponencialpdf(lambda, x);
Ecdf = exponencialcdf(lambda, x);

% Output:
% ExponencialPDF
subplot(2,1,1);
plot(x, Epdf, 'LineWidth', 2);
title('Exponencial PDF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
%xticks(x);
grid on;

% ExponencialCDF
subplot(2,1,2);
plot(x, Ecdf, 'LineWidth', 2);
title('Exponencial CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
%xticks(x);
grid on;
