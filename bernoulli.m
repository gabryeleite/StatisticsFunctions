clear
clc

% Funções:
function pv = bernoullipmf(p, x)
  % For Bernoulli (p) rv X
  % input = vector x
  % output = vector pv
  % such that pv(i) = Prob(X = x(i))
  % Na função bernoulli só são aceitos os valores 0 e 1;
  % portanto, qualquer valor diferente terá a probabilidade 0
  pv = (1 - p) * (x == 0) + p * (x == 1);
  pv = pv(:);
end

function cdf = bernoullicdf(pmf, x)
  % For Bernoulli (p) random variable X
  % input = vector x
  % output = vector cdf
  % such that cdf(i) = Prob(X <= x(i))
  x = floor(x(:));
  cdf = cumsum(pmf); % cumsum = soma acumulada das probabilidades
  cdf(x > 1) = 1; % valores de X maiores que 1 terão probabilidade 1
end

function x = bernoullirv(p, m)
  % return m samples of Bernoulli(p) random variable
  r = rand(m, 1);
  x = (r >= (1 - p));
endfunction

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

