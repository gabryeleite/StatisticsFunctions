clear
clc

% Funções:
% Bernoulli:
function pv = bernoullipmf(p, x)
  % Na função bernoulli só são aceitos os valores 0 e 1;
  % portanto, qualquer valor diferente terá a probabilidade 0
  pv = (1 - p) * (x == 0) + p * (x == 1);
  pv = pv(:);
end

function cdf = bernoullicdf(pmf, x)
  % soma cumulativa
  x = floor(x(:));
  cdf = cumsum(pmf); % cumsum = soma acumulada das probabilidades
  cdf(x > 1) = 1; % valores de X maiores que 1 terão probabilidade 1
end

function x = bernoullirv(p, m)
  % retornar m amostras da variável aleatória Bernoulli (p)
  r = rand(m, 1);
  x = (r >= (1 - p));
endfunction

function bernoulli()
  % Input:
  p = input('Digite a probabilidade p (valores entre 0 e 1): ');
  %m = input('Digite o numero de amostras: ');
  fid = fopen('dados.txt','r');
  x = fscanf(fid, '%f');
  fclose(fid);
  %x = bernoullirv(p, m);
  Bpmf = bernoullipmf(p, x); % bernoullipmf
  Bcdf = bernoullicdf(Bpmf, x); % bernoullicdf

  % Output:
  % BernoulliPMF
  subplot(2,1,1);
  stem(x, Bpmf, 'LineWidth', 2);
  title('Bernoulli PMF');
  xlabel('X');
  ylabel('Probabilidade');
  ylim([0, 1]);
  xticks(x);
  grid on;

  % BernoulliCDF
  subplot(2,1,2);
  stairs(x, Bcdf, 'LineWidth', 2); % stairs melhor representação para CDF, como se trata de uma soma acumulada
  title('Bernoulli CDF');
  xlabel('X');
  ylabel('Probabilidade acumulada');
  ylim([0, 1.5]);
  xticks(x);
  grid on;
end

% Binomial:
function pmf = binomialpmf(n, p, x)
  % x = vetor de inteiros não negativos
  if p < 0.5
    pp = p;
  else
    pp = 1 - p;
  end

  i = 0:n-1;
  ip = ((n - i) ./ (i + 1)) * (pp / (1 - pp));
  pb = ((1 - pp)^n) * cumprod([1, ip]);

  if pp < p
    pb = fliplr(pb);
  end

  pb = pb(:); % pb = [P[X=0] ... P[X=n]]^t
  x = x(:);
  okx = (x >= 0) .* (x <= n) .* (x == floor(x));
  x = okx .* x;
  pmf = okx .* pb(x + 1);
end

function cdf = binomialcdf(pmf, x)
  x = floor(x(:)); % para x(i) não inteiros
  allx = 0:max(x);

  % Calcula cdf de 0 até max(x)
  allcdf = cumsum(pmf);

  okx = (x >= 0); % x(i) < 0 probabilidade zero
  x = (okx .* x); % define probabilidade zero x(i) = 0

  cdf = okx .* allcdf(x + 1); % zero para zero-prob x(i)
end

function result = count(cdf, values)
  % Função de contagem para encontrar o primeiro índice onde cdf é maior que cada valor

  result = zeros(size(values));
  for i = 1:numel(values)
    result(i) = find(cdf >= values(i), 1, 'first') - 1;
  end
end

function x = binomialrv(n, p, m)
  % Gera m amostras da distribuição binomial(n, p)

  r = rand(m, 1);
  cdf = binomialcdf(n, p, 0:n);
  x = count(cdf, r);
end

function binomial()
  % Input:
  p = input('Digite a probabilidade p (valores entre 0 e 1): ');
  n = input('Digite o numero de ensaios n: ');
  %m = input('Digite o numero de amostras: ');
  fid = fopen('dados.txt','r');
  x = fscanf(fid, '%d');
  fclose(fid);
  %x = binomialrv(n, p, m);
  Bipmf = binomialpmf(n, p, x); % binomialpmf
  Bicdf = binomialcdf(Bipmf, x); % binomialcdf

  % Output:
  % BinomialPMF
  subplot(2,1,1);
  stem(x, Bipmf, 'LineWidth', 2);
  title('Binomial PMF');
  xlabel('X');
  ylabel('Probabilidade');
  ylim([0, 1]);
  xticks(x);
  grid on;

  % BinomialCDF
  subplot(2,1,2);
  stairs(x, Bicdf, 'LineWidth', 2); % stairs melhor representação para CDF, como se trata de uma soma acumulada
  title('Binomial CDF');
  xlabel('X');
  ylabel('Probabilidade acumulada');
  ylim([0, 1.5]);
  xticks(x);
  grid on;
end

% Exponencial:
function cdf = exponencialcdf(lambda, x)
    cdf = 1.0 - exp(-lambda * x);
end

function pdf = exponencialpdf(lambda, x)
    pdf = lambda * exp(-lambda * x) .* (x >= 0);
end

function rv = exponencialrv(lambda, m)
    rv = -(1/lambda) * log(1 - rand(m, 1));
end

function exponencial()
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
end

% Finita:
function pmf = finitepmf(sx, px, x)
    % pmf(i) = P[X = x(i)]

    pmf = zeros(size(x(:)));
    for i = 1:length(x)
        pmf(i) = sum(px(sx == x(i)));
    end
end

function cdf = finitecdf(s, p, x)
    % cdf(i) = P[X <= x(i)]

    cdf = zeros(size(x));
    for i = 1:length(x)
        pxi = sum(p(s <= x(i)));
        cdf(i) = pxi;
    end
end

function rho = finitecoeff(SX, SY, PXY)
    % Calcula o coeficiente de correlação rho de variáveis aleatórias finitas X e Y
    ex = finiteexp(SX, PXY);
    vx = finitevar(SX, PXY);
    ey = finiteexp(SY, PXY);
    vy = finitevar(SY, PXY);
    R = finiteexp(SX .* SY, PXY);
    rho = (R - ex * ey) / sqrt(vx * vy);
end

function covxy = finitecov(SX, SY, PXY)
    % Retorna a covariância das variáveis aleatórias finitas X e Y dadas pelas grades SX, SY e PXY
    ex = finiteexp(SX, PXY);
    ey = finiteexp(SY, PXY);
    R = finiteexp(SX .* SY, PXY);
    covxy = R - ex * ey;
end

function ex = finiteexp(sx, px)
    % Retorna o valor esperado E[X]
    ex = sum((sx(:)) .* (px(:)));
end

function v = finitevar(sx, px)
    % Retorna a variância Var[X]
    ex2 = finiteexp(sx.^2, px);
    ex = finiteexp(sx, px);
    v = ex2 - (ex^2);
end

function x = finiterv(s, p, m)
    % Retorna m amostras de uma variável aleatória finita (s, p)
    r = rand(m, 1);
    x = zeros(m, 1);

    for i = 1:m
        for j = 1:length(s)
            if r(i) <= sum(p(1:j))
                x(i) = s(j);
                break;
            end
        end
    end
end

function finita()
  %Input:
  % Definindo a grade de amostras e probabilidades
  %m = input('Digite o numero de amostras: '); % Número de amostras
  % sx = [0, 1, 2, 3]; % Exemplo de espaço amostral
  fid = fopen('dados.txt','r');
  sx = fscanf(fid, '%d'); % lendo o espaço amostral de um arquivo
  fclose(fid);
  % px = [0.1, 0.2, 0.3, 0.4]; % Exemplo de probabilidades
  fid = fopen('probabilidades.txt','r'); % lendo as probabilidades de um arquivo
  px = fscanf(fid, '%f'); % a soma das probabilidades tem que dar 1
  fclose(fid);            % e a quantidade deve ser a mesma do espaço amostral
  %x_values = finiterv(sx, px, m); % Gerando amostras da variável aleatória finita
  pmf = finitepmf(sx, px, sx); % Usei sx como valores de x para a PMF
  cdf = finitecdf(sx, px, sx); % Usei sx como valores de x para a CDF

  % Output:
  % FinitaPMF
  subplot(2, 1, 1);
  stem(sx, pmf, 'LineWidth', 2);
  title('Finita PMF');
  xlabel('X');
  ylabel('Probabilidade');
  ylim([0, 1]);
  xticks(sx);
  grid on;

  % FinitaCDF
  subplot(2, 1, 2);
  stairs(sx, cdf, 'LineWidth', 2);
  title('Finita CDF');
  xlabel('X');
  ylabel('Probabilidade acumulada');
  ylim([0, 1.5]);
  xticks(sx);
  grid on;
end

% Pascal:
function pmf = pascalpmf(k, p, x)
  % vetor pmf: pmf(i) = Prob[X = x(i)]

  x = x(:);
  n = max(x);
  i = (k:n-1)';
  ip = [1; (1-p) * (i./(i+1-k))];

  % pb = all n-k+1
  pb = (p^k) * cumprod(ip);

  okx = (x == floor(x)) & (x >= k);

  x = (okx .* x) + k * (1 - okx);

  % pmf(i) = 0 a menos que x(i) >= k
  pmf = okx .* pb(x-k+1);
end

function cdf = pascalcdf(k, p, x)
  % cdf(i) = Prob[X <= x(i)]

  x = floor(x(:));
  allx = k:max(x);

  allcdf = cumsum(pascalpmf(k, p, allx));

  okx = (x >= k);

  x = (okx .* x) + ((1 - okx) * k);
  cdf = okx .* allcdf(x - k + 1);
end

function x = pascalrv(k, p, m)
  % returna m amostras de pascal(k,p) rv

  r = rand(m, 1);
  rmax = max(r);
  xmin = k;
  xmax = ceil(2 * (k / p));

  sx = xmin:xmax;
  cdf = pascalcdf(k, p, sx);

  while cdf(end) <= rmax
      xmax = 2 * xmax;
      sx = xmin:xmax;
      cdf = pascalcdf(k, p, sx);
  end

  x = xmin + sum(cdf <= r, 2);
end

function pascal()
  % Input:
  k = input('Informe o parâmetro k (inteiro positivo): ');
  p = input('Informe a probabilidade p (entre 0 e 1): ');
  % m = input('Informe o número de amostras (m): ');
  fid = fopen('dados.txt','r');
  x = fscanf(fid, '%d');
  fclose(fid);
  %x = pascalrv(k, p, m)
  valores_pmf = pascalpmf(k, p, x);
  valores_cdf = pascalcdf(k, p, x);

  % Output:
  % PascalPMF
  subplot(2, 1, 1);
  stem(x, valores_pmf, 'LineWidth', 2);
  title('Pascal PMF');
  xlabel('X');
  ylabel('Probabilidade');
  ylim([0, 1]);
  xticks(x);
  grid on;

  % PascalCDF
  subplot(2, 1, 2);
  stairs(x, valores_cdf, 'LineWidth', 2);
  title('Pascal CDF');
  xlabel('X');
  ylabel('Probabilidade acumulada');
  ylim([0, 1.5]);
  xticks(x);
  grid on;
end

% Poisson:
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

function poisson()
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
end

% Programa principal:
while true
  clc
  % Exibir opções do menu
  fprintf('[1] Bernoulli\n');
  fprintf('[2] Binomial\n');
  fprintf('[3] Uniforme\n');
  fprintf('[4] Earlang\n');
  fprintf('[5] Exponencial\n');
  fprintf('[6] Finita\n');
  fprintf('[7] Gauss\n');
  fprintf('[8] Geometrica\n');
  fprintf('[9] Pascal\n');
  fprintf('[10] Poisson\n');
  fprintf('[0] Sair\n');
  opcao = input('Escolha: ');

  switch opcao
      case 1
          bernoulli();
      case 2
          binomial();
      case 5
          exponencial();
      case 6
          finita();
      case 9
          pascal();
      case 10
          poisson();
      case 0
          fprintf('Programa finalizado!\n');
          break; % Interrompe o loop
      otherwise
          fprintf('Opcao invalida! Tente novamente\n');
  end
end
