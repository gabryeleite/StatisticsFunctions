clear
clc

function pmf = duniformpmf(k, l, x)
  % discrete uniform(k,l) rv X,
  % input = vector x
  % output = vector pmf: pmf(i) = Prob[X = x(i)]

  pmf = (x >= k) .* (x <= l) .* (x == floor(x));
  pmf = pmf(:) / (l - k + 1);
end

function cdf = duniformcdf(pmf, x)
  % cdf(i) = Prob[X <= x(i)]
  x = floor(x(:)); 

  allcdf = cumsum(pmf);

  okx = (x >= 0); % x(i) < 0 probabilidade zero
  x = (okx .* x); % define probabilidade zero x(i) = 0

  cdf = okx .* allcdf(x + 1); % zero para zero-prob x(i)
end

function x = duniformrv(k, l, m)
  % returns m amostras
  r = rand(m, 1);
  cdf = duniformcdf(k, l, k:l);
  [~, x] = histc(r, [0; cumsum(cdf)]);
  x = x + k - 1;
end

% Input:
fprintf('Os limites devem estar de acordo com os dados passados para funcionar devidamente\n');
k = input('Digite o limite inferior: '); % Limite inferior da distribuição
l = input('Digite o limite superior: '); % Limite superior da distribuição
%m = input('Digite o numero de amostras: ');
%x = duniformrv(k, l, m);
fid = fopen('dados.txt','r');
x = fscanf(fid, '%d');
fclose(fid);
% Calcular PMF e CDF
valores_pmf = duniformpmf(k, l, x);
valores_cdf = duniformcdf(valores_pmf, x);

% UniformePMF
subplot(2, 1, 1);
stem(x, valores_pmf, 'LineWidth', 2);
title('Uniforme PMF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
xticks(x);
grid on;

% UniformeCDF
subplot(2, 1, 2);
stairs(x, valores_cdf, 'LineWidth', 2);
title('Uniforme CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
xticks(x);
grid on;
