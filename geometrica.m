clear
clc

% Funções:
function pmf = geometricapmf(p, x)
  % pmf(i) = Prob[X=x(i)]

  x = x(:);
  pmf = p * ((1-p).^(x-1));
  pmf = (x > 0) .* (x == floor(x)) .* pmf;
end

function cdf = geometricacdf(p, x)
  % for geometric(p) rv X,
  % For input vector x, output is vector
  % cdf such that cdf_i = Prob(X <= x_i)

  x = (x(:) >= 1) .* floor(x(:));
  cdf = 1 - ((1 - p).^x);
end

function x = geometricarv(p, m)
  % Usage: x = geometricrv(p, m)
  % Returns m samples of a geometric (p) rv

  r = rand(m, 1);
  x = ceil(log(1-r) / log(1-p));
end

% Input:
p = input('Digite a probabilidade p (entre 0 e 1): ');
%m = input('Digite o numero de amostras: ');
fid = fopen('dados.txt','r');
x = fscanf(fid, '%d');
fclose(fid);
%x = geometricarv(p, m);
Gepmf = geometricapmf(p, x); % geometricpmf
Gecdf = geometricacdf(p, x); % geometriccdf

% Output:
% GeometricaPMF
subplot(2,1,1);
stem(x, Gepmf, 'LineWidth', 2);
title('Geometrica PMF');
xlabel('X');
ylabel('Probabilidade');
ylim([0, 1]);
xticks(x);
grid on;

% GeometricaCDF
subplot(2,1,2);
stairs(x, Gecdf, 'LineWidth', 2); % stairs melhor representação para CDF, como se trata de uma soma acumulada
title('Geometrica CDF');
xlabel('X');
ylabel('Probabilidade acumulada');
ylim([0, 1.5]);
xticks(x);
grid on;
