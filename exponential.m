clear
clc

% Funções:
function cdf = exponentialcdf(lambda, x)
    cdf = 1.0 - exp(-lambda * x);
end

function pdf = exponentialpdf(lambda, x)
    pdf = lambda * exp(-lambda * x) .* (x >= 0);
end

function rv = exponentialrv(lambda, m)
    rv = -(1/lambda) * log(1 - rand(m, 1));
end

% Input:
lambda = input('Digite a taxa de chegada lambda (valor positivo): ');
m = input('Digite o numero de amostras: ');

% Gerando amostras exponenciais:
Erv = exponentialrv(lambda, m);
x_values = linspace(0, max(Erv), 1000);

% Output:
% ExponentialPDF
subplot(2,1,1);
Epdf = exponentialpdf(lambda, x_values);
plot(x_values, Epdf);
title('Exponential PDF');
xlabel('X');
ylabel('Densidade de Probabilidade');
grid on;

% ExponentialCDF
subplot(2,1,2);
Ecdf = exponentialcdf(lambda, x_values);
plot(x_values, Ecdf);
title('Exponential CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
grid on;
