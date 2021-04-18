clc;
clear;
close all;

rng(1);

n = 128;
tol=0.1;
rmse_avg = [];

for alpha=[0, 3]
    U = zeros(n);
    vi = randn(n,1);
    U(:,1) = vi ./ norm(vi);
    for i=2:n
      nrm = 0;
      while nrm < tol
        vi = randn(n,1);
        vi = vi - U(:,1:i-1) * (U(:,1:i-1)' * vi);
        nrm = norm(vi);
      end
      U(:,i) = vi ./ nrm;
    end

    Sigma = U * diag((1:n).^(-alpha)) * U';

    xs = mvnrnd(zeros(n,1), Sigma, 10)';

    for m=[40, 50, 64, 80, 100, 120]
        rmse = zeros(10, 1);
        phi = randn(m, n) / sqrt(m);

        for i=1:10
            y = phi * xs(:, i);
            sig = 0.01 * mean(abs(y));
            y = y + mvnrnd(zeros(m, 1), eye(m)*sig^2, 1)';
            x = (Sigma - Sigma*phi'/(eye(m)*sig^2 + phi*Sigma*phi')*phi*Sigma) * phi' * y / (sig^2);

            rmse(i) = sqrt(mean((x - xs(:, i)).^2)) / sqrt(mean(xs(:, i).^2));
        end
        rmse_avg = [rmse_avg mean(rmse)];
    end
end

