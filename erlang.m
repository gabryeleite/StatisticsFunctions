clear
clc

function f = erlangpmf(n, lambda, x)
    f = ((lambda^n) / factorial(n)) .* (x.^(n-1)) .* exp(-lambda*x);
end

function F = erlangcdf(n, lambda, x)
    F = 1 - gammainc(x * lambda, n, 'upper');
end

function pb = erlangb(rho, c)
  pn = exp(-rho) * poisspdf(0:c, rho);
  pb = pn(c + 1) / sum(pn);
end

function x = erlangrv(n, lambda, m)
    y = exponentialrv(lambda, m * n);
    y_reshaped = reshape(y, m, n);
    x = sum(y_reshaped, 2);
end

% Input:
n = input('Informe o parâmetro n (inteiro positivo): ');
lambda = input('Informe lambda (entre 0 e 1): ');
% m = input('Informe o número de amostras (m): ');
fid = fopen('dados.txt','r');
x = fscanf(fid, '%d');
fclose(fid);
%x = erlangrv(n, lambda, m)
valores_pmf = erlangpmf(n, lambda, x);
valores_cdf = erlangcdf(n, lambda, x);

% Output:
% ErlangPMF
subplot(2, 1, 1);
stem(x, valores_pmf, 'LineWidth', 2);
title('Erlang PMF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
xticks(x);
grid on;

% ErlangCDF
subplot(2, 1, 2);
stairs(x, valores_cdf, 'LineWidth', 2);
title('Erlang CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
xticks(x);
grid on;
