clear
clc

function pmf = poisson_pmf(alpha, x)
  % Variável aleatória X seguindo a distribuição de Poisson(alpha)
  % Saída: vetor pmf: pmf(i) = P[X = x(i)]

  x = x(:);
  k = (1:max(x))';
  logfatoriais = cumsum(log(k));
  pb = exp([-alpha; -alpha + (k * log(alpha)) - logfatoriais]);

  okx = (x >= 0) .* (x == floor(x));
  x = okx .* x;
  pmf = okx .* pb(x + 1); % pmf(i) = 0 para x(i) com probabilidade zero
end

function cdf = poisson_cdf(pmf, x)
  % Saída cdf(i) = Prob[X <= x(i)]

  x = floor(x(:));
  sx = 0:max(x);

  allcdf = cumsum(pmf); % soma cumulativa da PMF para obter a CDF

  % Ajuste para garantir que a CDF atinja 1 no final
  %cdf = cdf / cdf(end);

  % Correção para x(i) < 0
  okx = (x >= 0);
  cdf = okx .* allcdf(x + 1); % cdf = 0 para x(i) < 0
end

function x = variavel_poisson(alpha, m)
  % Retorna m amostras da variável aleatória de Poisson(alpha) X
  r = rand(m, 1);
  rmax = max(r);
  xmin = 0;
  xmax = ceil(2 * alpha); % define a faixa máxima
  sx = xmin:xmax;
  cdf = poisson_cdf(alpha, sx);

  % enquanto sum(cdf <= rmax) == (xmax - xmin + 1)
  while cdf(end) <= rmax
    xmax = 2 * xmax;
    sx = xmin:xmax;
    cdf = poisson_cdf(alpha, sx);
  end

  x = xmin + sum(cdf <= rmax);
end

% Input:
alpha = input('Informe o parâmetro alpha para a distribuição de Poisson: ');
% m = input('Informe o número de amostras (m): ');
fid = fopen('dados.txt','r');
valores_x = fscanf(fid, '%d');
fclose(fid);
%valores_x = variavel_poisson(alpha, m) % Gera variáveis aleatórias seguindo a distribuição de Poisson
valores_pmf = poisson_pmf(alpha, valores_x);
valores_cdf = poisson_cdf(valores_pmf, valores_x);

% Output:
% PoissonPMF
% figure;
subplot(2, 1, 1);
stem(valores_x, valores_pmf, 'LineWidth', 2);
title('Poisson PMF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
xticks(valores_x);
grid on;

% PoissonCDF
subplot(2, 1, 2);
stairs(valores_x, valores_cdf, 'LineWidth', 2);
title('Poisson CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
xticks(valores_x);
grid on;
