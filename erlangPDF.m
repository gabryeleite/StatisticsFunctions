clear
clc

function f = erlangpdf(n, lambda, x)
    f = ((lambda^n) / factorial(n)) .* (x.^(n-1)) .* exp(-lambda*x);
end

function F = erlangcdf(n, lambda, x)
    F = 1 - gammainc(x * lambda, n, 'upper');
end

function x = erlangrv(n, lambda, m)
    y = exponentialrv(lambda, m*n);
    x = sum(reshape(y, m, n), 2);
end

% Input:
a = input('Digite o limite inferior: '); % Limite inferior da distribuição
b = input('Digite o limite superior: '); % Limite superior da distribuição
n = input('Informe o parâmetro n (inteiro positivo): ');
lambda = input('Informe lambda (entre 0 e 1): ');
% m = input('Informe o número de amostras (m): ');
%x = erlangrv(n, lambda, m);
x = a : 0.1 : b;

%fid = fopen('dados.txt','r');
%x = fscanf(fid, '%d');
%fclose(fid);

% Calcular PDF e CDF
valores_pdf = erlangpdf(n, lambda, x);
valores_cdf = erlangcdf(n, lambda, x);

% Output:
% ErlangPDF
subplot(2, 1, 1);
plot(x, valores_pdf, 'LineWidth', 2);
title('Erlang PDF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
%xticks(x);
grid on;

% ErlangCDF
subplot(2, 1, 2);
plot(x, valores_cdf, 'LineWidth', 2);
title('Erlang CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
%xticks(x);
grid on;
