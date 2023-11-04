clear
clc

% Funções:
function pmf = binomialpmf(n, p, x)
  % binomial(n, p) rv X,
  % input = vector x
  % output = vector pmf: pmf(i) = Prob[X = x(i)]

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
  % Usage: cdf = binomialcdf(n, p, x)
  % For binomial(n, p) rv X,
  % and input vector x, output is
  % vector cdf: cdf(i) = P[X <= x(i)]

  x = floor(x(:)); % for non-integer x(i)
  allx = 0:max(x);

  % Calculate cdf from 0 to max(x)
  allcdf = cumsum(pmf);

  okx = (x >= 0); % x(i) < 0 are zero-prob values
  x = (okx .* x); % set zero-prob x(i) = 0

  cdf = okx .* allcdf(x + 1); % zero for zero-prob x(i)
end

function result = count(cdf, values)
  % Count function to find the first index where cdf is greater than each value

  result = zeros(size(values));
  for i = 1:numel(values)
    result(i) = find(cdf >= values(i), 1, 'first') - 1;
  end
end

function x = binomialrv(n, p, m)
  % Usage: x = binomialrv(n, p, m)
  % Generate m samples from the binomial(n, p) distribution

  r = rand(m, 1);
  cdf = binomialcdf(n, p, 0:n);
  x = count(cdf, r);
end

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
