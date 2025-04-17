load('I.mat');

dat = Idata.I(1:100:end);

z1test(dat)


function [kmedian, c, kcorr, Mall, p_best, q_best, M_best, t] = z1test(x)
    if size(x,2) == 1
        x = x';
    end

    N = length(x);
    j = 1:N;
    t = 1:round(N/10);
    num_c = 100;

    Mall = zeros(num_c, length(t));
    kcorr = zeros(1, num_c);
    c = pi/5 + rand(1, num_c) * (3*pi/5);

    % Preallocate for best iteration
    best_idx = 1;
    max_corr = -inf;

    for its = 1:num_c
        p = cumsum(x .* cos(j * c(its)));
        q = cumsum(x .* sin(j * c(its)));

        M = zeros(1, length(t));

        for n = t
            M(n) = mean((p(n+1:N) - p(1:N-n)).^2 + ...
                        (q(n+1:N) - q(1:N-n)).^2) - ...
                   mean(x)^2 * (1 - cos(n * c(its))) / (1 - cos(c(its)));
        end

        Mall(its, :) = M;
        kcorr(its) = corr(t', M');

        if kcorr(its) > max_corr
            max_corr = kcorr(its);
            best_idx = its;
            p_best = p;
            q_best = q;
            M_best = M;
        end
    end

    if (max(x) - min(x)) / mean(abs(diff(x))) > 10 || ...
       median(kcorr(c < mean(c))) - median(kcorr(c > mean(c))) > 0.5
        disp('Warning: data is probably oversampled.')
        disp('Use coarser sampling or reduce the maximum value of c.')
    end

    kmedian = median(kcorr);
end

load('I.mat');

[kmedian, c, kcorr, Mall, p, q, M, t] = z1test(Idata.I);

% Plot p vs. q
figure;
plot(p, q, 'LineWidth',2);
xlabel('p');
ylabel('q');
title('Trajectory: q vs. p');

% Plot M vs. t
figure;
plot(t, M, 'LineWidth',2);
xlabel('t');
ylabel('M');
title('M vs. t');

figure;
subplot(2,1,1);
plot(Idata.I); title('Original Signal');

subplot(2,1,2);
plot(dat); title('Downsampled Signal');