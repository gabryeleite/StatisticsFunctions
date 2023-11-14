clear
clc

function pmf = pascalpmf(k, p, x)
    % For Pascal (k,p) rv X, and
    % input vector x, output is a
    % vector pmf: pmf(i) = Prob[X = x(i)]

    x = x(:);
    n = max(x);
    i = (k:n-1)';
    ip = [1; (1-p) * (i./(i+1-k))];

    % pb = all n-k+1 pascal probs
    pb = (p^k) * cumprod(ip);

    okx = (x == floor(x)) & (x >= k);
    % set bad x(i) = k to stop bad indexing
    x = (okx .* x) + k * (1 - okx);

    % pmf(i) = 0 unless x(i) >= k
    pmf = okx .* pb(x-k+1);
end

function cdf = pascalcdf(k, p, x)
    % Usage: cdf = pascalcdf(k, p, x)
    % For a Pascal (k, p) rv X
    % and input vector x, the output
    % is a vector cdf such that
    % cdf(i) = Prob[X <= x(i)]

    x = floor(x(:)); % for non-integer x(i)
    allx = k:max(x);
    % allcdf holds all needed cdf values
    allcdf = cumsum(pascalpmf(k, p, allx));

    % x_i < k have zero-prob,
    % other values are OK
    okx = (x >= k);

    % set zero-prob x(i) = k,
    % just so indexing is not fouled up
    x = (okx .* x) + ((1 - okx) * k);
    cdf = okx .* allcdf(x - k + 1);
end

function x = pascalrv(k, p, m)
    % return m samples of pascal(k,p) rv

    r = rand(m, 1);
    rmax = max(r);
    xmin = k;
    xmax = ceil(2 * (k / p)); % set max range

    sx = xmin:xmax;
    cdf = pascalcdf(k, p, sx);

    while cdf(end) <= rmax
        xmax = 2 * xmax;
        sx = xmin:xmax;
        cdf = pascalcdf(k, p, sx);
    end

    x = xmin + sum(cdf <= r, 2);
end

% Input:
k = input('Informe o parâmetro k (inteiro positivo): ');
p = input('Informe a probabilidade p (entre o e 1): ');
% m = input('Informe o número de amostras (m): ');
fid = fopen('dados.txt','r');
x = fscanf(fid, '%d');
fclose(fid);
%x = pascalrv(k, p, m) % Gera variáveis aleatórias seguindo a distribuição de Poisson
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
