clear
clc

% Funções:
function cdf = exponencialcdf(lambda, x)
    cdf = 1.0 - exp(-lambda * x);
end

function pdf = exponencialpdf(lambda, x)
    pdf = lambda * exp(-lambda * x) .* (x >= 0);
end

function rv = exponencialrv(lambda, m)
    rv = -(1/lambda) * log(1 - rand(m, 1));
end

% Input:
lambda = input('Digite a taxa de chegada lambda (entre 0 e 1): ');
%m = input('Digite o numero de amostras: ');
fid = fopen('dados.txt','r');
x = fscanf(fid, '%d');
fclose(fid);
% Gerando amostras exponenciais:
%Erv = exponencialrv(lambda, m);
%x = linspace(0, max(Erv), m); % 1000
Epdf = exponencialpdf(lambda, x);
Ecdf = exponencialcdf(lambda, x);

% Output:
% ExponencialPDF
subplot(2,1,1);
stem(x, Epdf, 'LineWidth', 2);
title('Exponencial PDF');
xlabel('X');
ylabel('Densidade de Probabilidade');
ylim([0, 1]);
xticks(x);
grid on;

% ExponencialCDF
subplot(2,1,2);
stairs(x, Ecdf, 'LineWidth', 2);
title('Exponencial CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
xticks(x);
grid on;
