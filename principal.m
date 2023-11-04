clear
clc

% Funções:
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
  % Gera m amostras da variável aleatória Bernoulli (p)
  r = rand(m, 1);
  x = (r >= (1 - p));
endfunction

function bernoulli()
  % Input:
  p = input('Digite a probabilidade p (valores entre 0 e 1): ');
  %m = input('Digite o numero de amostras: ');
  fid = fopen('bernoulli.txt','r');
  x = fscanf(fid, '%f');
  fclose(fid);
  %x = bernoullirv(p, m);
  Bpmf = bernoullipmf(p, x); % bernoullipmf
  Bcdf = bernoullicdf(Bpmf, x); % bernoullicdf

  % Output:
  % BernoulliPMF
  subplot(2,1,1);
  stem(x, Bpmf);
  title('Bernoulli PMF');
  xlabel('X');
  ylabel('Probabilidade');
  ylim([0, 1]);
  xticks(x);
  grid on;

  % BernoulliCDF
  subplot(2,1,2);
  stairs(x, Bcdf); % stairs melhor representação para CDF, como se trata de uma soma acumulada
  title('Bernoulli CDF');
  xlabel('X');
  ylabel('Probabilidade acumulada');
  ylim([0, 1.5]);
  xticks(x);
  grid on;
end

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
  fid = fopen('binomial.txt','r');
  x = fscanf(fid, '%d');
  fclose(fid);
  %x = binomialrv(n, p, m);
  Bipmf = binomialpmf(n, p, x); % binomialpmf
  Bicdf = binomialcdf(Bipmf, x); % binomialcdf

  % Output:
  % BinomialPMF
  subplot(2,1,1);
  stem(x, Bipmf);
  title('Binomial PMF');
  xlabel('X');
  ylabel('Probabilidade');
  ylim([0, 1]);
  xticks(x);
  grid on;

  % BinomialCDF
  subplot(2,1,2);
  stairs(x, Bicdf); % stairs melhor representação para CDF, como se trata de uma soma acumulada
  title('Binomial CDF');
  xlabel('X');
  ylabel('Probabilidade acumulada');
  ylim([0, 1.5]);
  xticks(x);
  grid on;
end

% Programa principal:
while true
      % Exibir opções do menu
      fprintf('[1] Bernoulli\n');
      fprintf('[2] Binomial\n');
      fprintf('[0] Sair\n');
      opcao = input('Escolha: ');

      switch opcao
          case 1
              bernoulli();
          case 2
              binomial();
          case 0
              fprintf('Programa finalizado!\n');
              break; % Interrompe o loop
          otherwise
              fprintf('Opcao invalida! Tente novamente\n');
      end
end
