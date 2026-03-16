%% Chi-Square Goodness-of-Fit Tutorial
% This script walks through the mechanics of chi-square (chi^2)
% goodness-of-fit testing: what it measures, how to compute it, and how to
% interpret the result using the reduced chi-square and p-value.
%
% Key idea: if your model is correct and your error bars are honest, the
% chi^2 statistic follows a known distribution (the chi-squared
% distribution with nu = N - M degrees of freedom, where N is the number
% of data points and M is the number of fitted parameters).

%% 1. Define the true model
% We use a damped sinusoid:  y(x) = A*sin(k*x + phi)*exp(-gamma*x)
% This has four parameters: amplitude A, wavenumber k, phase phi, and
% damping rate gamma.

clear; close all; rng(42);  % reproducibility

% --- "Hidden" true parameters (what we will try to recover) ---
A_true     = 3.0;
k_true     = 2.5;
phi_true   = 0.8;
gamma_true = 0.15;

p_true = [A_true, k_true, phi_true, gamma_true];

% Model function handle: p = [A, k, phi, gamma]
model = @(x, p) p(1) * sin(p(2)*x + p(3)) .* exp(-p(4)*x);

%% 2. Generate synthetic data with Gaussian noise
N      = 60;                        % number of data points
x_data = linspace(0, 10, N)';      % evenly spaced x values
sigma  = 0.5;                       % known measurement uncertainty (same for all points)

y_true = model(x_data, p_true);
noise  = sigma * randn(N, 1);      % Gaussian noise ~ N(0, sigma)
y_data = y_true + noise;

figure(1); clf;
errorbar(x_data, y_data, sigma*ones(N,1), 'o', 'MarkerSize', 4, ...
    'DisplayName', 'Data (with error bars)');
hold on;
x_fine = linspace(0, 10, 500)';
plot(x_fine, model(x_fine, p_true), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'True model');
xlabel('x'); ylabel('y');
title('Synthetic Data from Damped Sinusoid');
legend('Location', 'best');
grid on;

%% 3. Fit the model to the data using least-squares (minimizing chi^2)
% The chi-square statistic is:
%
%   chi^2 = sum_i [ (y_i - model(x_i, p))^2 / sigma_i^2 ]
%
% Minimizing chi^2 over the parameters p gives the best-fit parameters.
% We use MATLAB's lsqcurvefit, which minimizes sum of squared residuals.
% Since all sigma_i are equal here, minimizing sum((y-model)^2) is
% equivalent to minimizing chi^2.

% For lsqcurvefit we need: F(p,x) returns the model prediction
fitfun = @(p, x) p(1) * sin(p(2)*x + p(3)) .* exp(-p(4)*x);

% Initial guesses (deliberately offset from truth to show the fit works)
p0 = [2.0, 2.0, 0.0, 0.1];

% Lower and upper bounds (optional but helps convergence)
lb = [0,   0, -pi, 0];
ub = [10, 10,  pi, 2];

options = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4);

[p_fit, resnorm] = lsqcurvefit(fitfun, p0, x_data, y_data, lb, ub, options);

fprintf('\n=== Fitted vs. True Parameters ===\n');
param_names = {'A', 'k', 'phi', 'gamma'};
for i = 1:4
    fprintf('  %-6s  true = %7.4f   fit = %7.4f\n', ...
        param_names{i}, p_true(i), p_fit(i));
end

%% 4. Compute the chi-square statistic
% This is the central quantity for goodness-of-fit testing.

y_fit    = model(x_data, p_fit);
residual = y_data - y_fit;

chi2 = sum(residual.^2 / sigma^2);

fprintf('\n=== Chi-Square Statistic ===\n');
fprintf('  chi^2 = %.2f\n', chi2);

% Note: lsqcurvefit returns resnorm = sum(residual.^2), so
% chi2 should equal resnorm / sigma^2.  Let's verify:
fprintf('  Verification: resnorm/sigma^2 = %.2f\n', resnorm / sigma^2);

%% 5. Degrees of freedom and reduced chi-square
% The number of degrees of freedom is:
%   nu = N - M
% where N = number of data points, M = number of fitted parameters.
%
% The REDUCED chi-square is:
%   chi^2_red = chi^2 / nu
%
% Interpretation:
%   chi^2_red ~ 1   => good fit (model explains data within errors)
%   chi^2_red >> 1   => bad fit (model doesn't describe data, or errors
%                       are underestimated)
%   chi^2_red << 1   => errors are overestimated, or the model is
%                       "overfitting" the noise

M  = length(p_fit);   % number of free parameters
nu = N - M;           % degrees of freedom

chi2_red = chi2 / nu;

fprintf('\n=== Reduced Chi-Square ===\n');
fprintf('  N (data points)       = %d\n', N);
fprintf('  M (free parameters)   = %d\n', M);
fprintf('  nu = N - M            = %d\n', nu);
fprintf('  chi^2_red = chi^2/nu  = %.4f\n', chi2_red);
fprintf('  Expected value if model is correct: chi^2_red ~ 1 +/- %.4f\n', ...
    sqrt(2/nu));

%% 6. P-value from the chi-squared distribution
% The p-value answers: "If the model is correct and the errors are right,
% what is the probability of getting a chi^2 value THIS LARGE or larger?"
%
%   p-value = 1 - chi2cdf(chi^2, nu)
%
% Conventionally:
%   p > 0.05  =>  no reason to reject the model
%   p < 0.05  =>  the fit is suspiciously bad (or errors are wrong)

p_value = 1 - chi2cdf(chi2, nu);

fprintf('\n=== P-Value ===\n');
fprintf('  p-value = %.4f\n', p_value);
if p_value > 0.05
    fprintf('  => No reason to reject the model (p > 0.05).\n');
else
    fprintf('  => Model may be inadequate or errors may be wrong (p < 0.05).\n');
end

%% 7. Visualize the chi-squared distribution and our result
figure(2); clf;

chi2_vals = linspace(0, 2*nu, 500);
chi2_pdf  = chi2pdf(chi2_vals, nu);

plot(chi2_vals, chi2_pdf, 'b-', 'LineWidth', 1.5);
hold on;
xline(chi2, 'r--', 'LineWidth', 2, 'Label', sprintf('\\chi^2 = %.1f', chi2));
xline(nu, 'k:', 'LineWidth', 1, 'Label', sprintf('\\nu = %d (expected)', nu));

xlabel('\chi^2');
ylabel('Probability density');
title(sprintf('\\chi^2 distribution (\\nu = %d) and observed \\chi^2', nu));
grid on;

% Shade the rejection region (right tail beyond our chi^2)
idx_shade = chi2_vals >= chi2;
area(chi2_vals(idx_shade), chi2_pdf(idx_shade), ...
    'FaceColor', [1 0.6 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('p-value = %.3f', p_value));
legend('Location', 'best');

%% 8. Plot the fit result
figure(3); clf;

subplot(2,1,1);
errorbar(x_data, y_data, sigma*ones(N,1), 'o', 'MarkerSize', 4, ...
    'DisplayName', 'Data');
hold on;
plot(x_fine, model(x_fine, p_fit), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Best fit');
plot(x_fine, model(x_fine, p_true), 'g--', 'LineWidth', 1, ...
    'DisplayName', 'True model');
xlabel('x'); ylabel('y');
title(sprintf('Best Fit (\\chi^2_{red} = %.3f, p = %.3f)', chi2_red, p_value));
legend('Location', 'best');
grid on;

subplot(2,1,2);
errorbar(x_data, residual, sigma*ones(N,1), 'o', 'MarkerSize', 4);
hold on;
yline(0, 'k-');
yline(sigma, 'r--', '+\sigma');
yline(-sigma, 'r--', '-\sigma');
xlabel('x'); ylabel('Residual (data - fit)');
title('Residuals');
grid on;

%% 9. Demonstration: what happens with a WRONG model?
% Let's deliberately fit a model that cannot describe the data (a simple
% exponential decay, ignoring the oscillation) and see how chi^2 responds.

fprintf('\n========================================\n');
fprintf('=== Wrong Model Demonstration ===\n');
fprintf('========================================\n');

wrong_model = @(p, x) p(1) * exp(-p(2)*x);

p0_wrong = [3.0, 0.1];
lb_wrong = [0, 0];
ub_wrong = [10, 2];

[p_wrong, resnorm_wrong] = lsqcurvefit(wrong_model, p0_wrong, ...
    x_data, y_data, lb_wrong, ub_wrong, options);

chi2_wrong    = resnorm_wrong / sigma^2;
M_wrong       = length(p_wrong);
nu_wrong      = N - M_wrong;
chi2_red_wrong = chi2_wrong / nu_wrong;
p_value_wrong  = 1 - chi2cdf(chi2_wrong, nu_wrong);

fprintf('  Wrong model: y = A*exp(-gamma*x)\n');
fprintf('  chi^2       = %.2f\n', chi2_wrong);
fprintf('  chi^2_red   = %.4f   (should be ~1 for a good fit)\n', chi2_red_wrong);
fprintf('  p-value     = %.2e\n', p_value_wrong);
fprintf('  => chi^2_red >> 1 and p ~ 0: the model is clearly wrong.\n');

figure(4); clf;

subplot(2,1,1);
errorbar(x_data, y_data, sigma*ones(N,1), 'o', 'MarkerSize', 4, ...
    'DisplayName', 'Data');
hold on;
plot(x_fine, wrong_model(p_wrong, x_fine), 'm-', 'LineWidth', 1.5, ...
    'DisplayName', 'Wrong model fit');
plot(x_fine, model(x_fine, p_true), 'g--', 'LineWidth', 1, ...
    'DisplayName', 'True model');
xlabel('x'); ylabel('y');
title(sprintf('Wrong Model (\\chi^2_{red} = %.1f, p = %.1e)', ...
    chi2_red_wrong, p_value_wrong));
legend('Location', 'best');
grid on;

subplot(2,1,2);
residual_wrong = y_data - wrong_model(p_wrong, x_data);
errorbar(x_data, residual_wrong, sigma*ones(N,1), 'o', 'MarkerSize', 4);
hold on;
yline(0, 'k-');
xlabel('x'); ylabel('Residual');
title('Residuals — Wrong Model (note the structure!)');
grid on;

%% 10. Demonstration: what happens with WRONG error bars?
% Even with the correct model, if sigma is wrong, chi^2 will be off.

fprintf('\n========================================\n');
fprintf('=== Wrong Error Bars Demonstration ===\n');
fprintf('========================================\n');

% Case A: sigma too small (underestimated errors)
sigma_small = 0.2;
chi2_small     = sum(residual.^2 / sigma_small^2);
chi2_red_small = chi2_small / nu;
p_value_small  = 1 - chi2cdf(chi2_small, nu);

fprintf('  Correct model, but sigma = %.2f (true sigma = %.2f):\n', ...
    sigma_small, sigma);
fprintf('    chi^2_red = %.3f,  p-value = %.2e\n', chi2_red_small, p_value_small);
fprintf('    => chi^2_red >> 1: errors look underestimated.\n\n');

% Case B: sigma too large (overestimated errors)
sigma_large = 1.5;
chi2_large     = sum(residual.^2 / sigma_large^2);
chi2_red_large = chi2_large / nu;
p_value_large  = 1 - chi2cdf(chi2_large, nu);

fprintf('  Correct model, but sigma = %.2f (true sigma = %.2f):\n', ...
    sigma_large, sigma);
fprintf('    chi^2_red = %.3f,  p-value = %.4f\n', chi2_red_large, p_value_large);
fprintf('    => chi^2_red << 1: errors look overestimated.\n');

% --- Figure 5: Underestimated errors (sigma too small) ---
figure(5); clf;

subplot(2,1,1);
errorbar(x_data, y_data, sigma_small*ones(N,1), 'o', 'MarkerSize', 4, ...
    'DisplayName', 'Data (tiny error bars)');
hold on;
plot(x_fine, model(x_fine, p_fit), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Best fit (correct model)');
xlabel('x'); ylabel('y');
title(sprintf('Underestimated Errors (\\sigma = %.2f, true = %.2f)  \\chi^2_{red} = %.2f,  p = %.1e', ...
    sigma_small, sigma, chi2_red_small, p_value_small));
legend('Location', 'best');
grid on;

subplot(2,1,2);
errorbar(x_data, residual, sigma_small*ones(N,1), 'o', 'MarkerSize', 4);
hold on;
yline(0, 'k-');
yline(sigma_small, 'r--', ['+\sigma = ' num2str(sigma_small)]);
yline(-sigma_small, 'r--', ['-\sigma = ' num2str(sigma_small)]);
xlabel('x'); ylabel('Residual (data - fit)');
title('Residuals — many points fall outside error bars');
grid on;

% --- Figure 6: Overestimated errors (sigma too large) ---
figure(6); clf;

subplot(2,1,1);
errorbar(x_data, y_data, sigma_large*ones(N,1), 'o', 'MarkerSize', 4, ...
    'DisplayName', 'Data (huge error bars)');
hold on;
plot(x_fine, model(x_fine, p_fit), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Best fit (correct model)');
xlabel('x'); ylabel('y');
title(sprintf('Overestimated Errors (\\sigma = %.2f, true = %.2f)  \\chi^2_{red} = %.3f,  p = %.3f', ...
    sigma_large, sigma, chi2_red_large, p_value_large));
legend('Location', 'best');
grid on;

subplot(2,1,2);
errorbar(x_data, residual, sigma_large*ones(N,1), 'o', 'MarkerSize', 4);
hold on;
yline(0, 'k-');
yline(sigma_large, 'r--', ['+\sigma = ' num2str(sigma_large)]);
yline(-sigma_large, 'r--', ['-\sigma = ' num2str(sigma_large)]);
xlabel('x'); ylabel('Residual (data - fit)');
title('Residuals — all points well within error bars');
grid on;

%% 11. Profile chi-square: rigorous confidence intervals for a parameter
% The key idea: to find the confidence interval for parameter p_j, we
% fix p_j at a series of trial values and RE-FIT all other parameters.
% At each trial value we record the best achievable chi^2.  This gives
% us a "profile" chi^2(p_j).
%
% The minimum of the profile is the global best-fit chi^2_min (achieved
% at p_j = p_j_best).  The confidence interval is the region where:
%
%   Delta_chi^2 = chi^2(p_j) - chi^2_min  <=  threshold
%
% For a SINGLE parameter at 95% confidence, the threshold is:
%
%   Delta_chi^2_95 = chi2inv(0.95, 1) = 3.84
%
% This works because fixing one parameter removes one degree of freedom,
% so the change in chi^2 follows a chi-squared distribution with 1 dof
% (Wilks' theorem).  The 95% interval is where the fit is "not
% significantly worse" than the best fit at the 95% level.

fprintf('\n================================================\n');
fprintf('=== Profile Chi-Square Confidence Intervals ===\n');
fprintf('================================================\n');

% We'll profile the wavenumber k (parameter 2) as a demonstration.
% The same procedure applies to any parameter.

delta_chi2_95 = chi2inv(0.95, 1);  % = 3.84 for 1 parameter
fprintf('  Delta chi^2 threshold (95%%, 1 dof) = %.4f\n\n', delta_chi2_95);

% Set up a grid of trial k values around the best fit
k_best  = p_fit(2);
k_range = linspace(k_best - 0.4, k_best + 0.4, 80);

% For each trial k: fix k, fit the remaining 3 parameters (A, phi, gamma)
chi2_profile = zeros(size(k_range));

% Model with k fixed: free parameters are q = [A, phi, gamma]
model_k_fixed = @(q, x, k_fixed) q(1) * sin(k_fixed*x + q(2)) .* exp(-q(3)*x);

% Initial guesses and bounds for the 3 free parameters
q0 = [p_fit(1), p_fit(3), p_fit(4)];   % start from best-fit values
lb_q = [0, -pi, 0];
ub_q = [10, pi, 2];

for i = 1:length(k_range)
    k_trial = k_range(i);

    % Wrap the model so lsqcurvefit only sees the free parameters
    fitfun_fixed_k = @(q, x) model_k_fixed(q, x, k_trial);

    [~, resnorm_i] = lsqcurvefit(fitfun_fixed_k, q0, x_data, y_data, ...
        lb_q, ub_q, options);

    chi2_profile(i) = resnorm_i / sigma^2;
end

% The profile Delta chi^2
delta_chi2_profile = chi2_profile - chi2;   % chi2 is the global minimum

% Find the 95% confidence bounds by interpolation
% Left bound: find where delta_chi2 crosses the threshold on the left
idx_left  = find(k_range < k_best & delta_chi2_profile > delta_chi2_95, 1, 'last');
idx_right = find(k_range > k_best & delta_chi2_profile > delta_chi2_95, 1, 'first');

if ~isempty(idx_left)
    % Interpolate between the last point above and first point below
    k_lo = interp1(delta_chi2_profile(idx_left:idx_left+1), ...
                   k_range(idx_left:idx_left+1), delta_chi2_95);
else
    k_lo = k_range(1);
    fprintf('  Warning: left 95%% bound not captured; widen k_range.\n');
end

if ~isempty(idx_right)
    k_hi = interp1(delta_chi2_profile(idx_right-1:idx_right), ...
                   k_range(idx_right-1:idx_right), delta_chi2_95);
else
    k_hi = k_range(end);
    fprintf('  Warning: right 95%% bound not captured; widen k_range.\n');
end

fprintf('  Profiled parameter: k (wavenumber)\n');
fprintf('  Best-fit value:     k = %.4f\n', k_best);
fprintf('  True value:         k = %.4f\n', k_true);
fprintf('  95%% confidence interval: k in [%.4f, %.4f]\n', k_lo, k_hi);
fprintf('  Asymmetric error bars:   k = %.4f  +%.4f / -%.4f\n', ...
    k_best, k_hi - k_best, k_best - k_lo);

if k_true >= k_lo && k_true <= k_hi
    fprintf('  => True value IS inside the 95%% interval (as expected).\n');
else
    fprintf('  => True value is OUTSIDE the 95%% interval.\n');
end

% --- Figure 7: Profile chi-square curve ---
figure(7); clf;

subplot(2,1,1);
plot(k_range, chi2_profile, 'b-', 'LineWidth', 1.5);
hold on;
yline(chi2, 'k:', 'LineWidth', 1, 'Label', '\chi^2_{min}');
yline(chi2 + delta_chi2_95, 'r--', 'LineWidth', 1.5, ...
    'Label', '\chi^2_{min} + 3.84');
xline(k_best, 'k:', 'LineWidth', 1);
xline(k_lo, 'r-', 'LineWidth', 1, 'Label', 'k_{lo}');
xline(k_hi, 'r-', 'LineWidth', 1, 'Label', 'k_{hi}');
xline(k_true, 'g--', 'LineWidth', 1, 'Label', 'k_{true}');
xlabel('k (fixed)');
ylabel('\chi^2');
title('Profile \chi^2: fix k, re-fit all other parameters');
grid on;

subplot(2,1,2);
plot(k_range, delta_chi2_profile, 'b-', 'LineWidth', 1.5);
hold on;
yline(delta_chi2_95, 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('\\Delta\\chi^2 = %.2f  (95%%)', delta_chi2_95));
yline(0, 'k:');

% Shade the 95% confidence region
idx_in = delta_chi2_profile <= delta_chi2_95;
area(k_range(idx_in), delta_chi2_profile(idx_in), ...
    'FaceColor', [0.7 0.85 1.0], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

xline(k_best, 'k:', 'LineWidth', 1);
xline(k_true, 'g--', 'LineWidth', 1, 'Label', 'k_{true}');
xlabel('k (fixed)');
ylabel('\Delta\chi^2 = \chi^2(k) - \chi^2_{min}');
title(sprintf('95%% confidence interval: k \\in [%.4f, %.4f]', k_lo, k_hi));
grid on;

%% Summary
fprintf('\n========================================\n');
fprintf('=== Summary ===\n');
fprintf('========================================\n');
fprintf(['The chi-square statistic measures how well a model describes\n' ...
    'data, relative to the measurement uncertainties:\n\n' ...
    '  chi^2 = sum( (data_i - model_i)^2 / sigma_i^2 )\n\n' ...
    'Key diagnostics:\n' ...
    '  1. Reduced chi^2 = chi^2 / (N - M) should be ~1.\n' ...
    '  2. p-value = 1 - chi2cdf(chi^2, N-M) should be > 0.05.\n' ...
    '  3. chi^2_red >> 1 => bad model or underestimated errors.\n' ...
    '  4. chi^2_red << 1 => overestimated errors (or overfitting).\n' ...
    '  5. ALWAYS inspect residuals for systematic structure.\n' ...
    '  6. Profile chi^2 gives rigorous confidence intervals:\n' ...
    '     fix one parameter, re-fit the rest, and find where\n' ...
    '     Delta_chi^2 = chi^2(p_j) - chi^2_min exceeds the\n' ...
    '     threshold chi2inv(CL, 1)  (3.84 for 95%%, 1 param).\n']);
