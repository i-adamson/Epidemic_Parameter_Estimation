clc; clear; close all;

% Define the parameter space for beta and gamma
beta_values = linspace(0, 3, 100);  % Transmission rate (beta)
gamma_values = linspace(0, 3, 100);  % Recovery rate (gamma)

% Create a grid for beta and gamma
[Beta, Gamma] = meshgrid(beta_values, gamma_values);

% Preallocate I_star for storing the steady-state infected population
I_star = zeros(size(Beta));

% Total population (normalized to 1)
N = 1;
kappa = 1;  % Incubation rate (assumed constant)

% Compute I* for each (beta, gamma) pair
for i = 1:length(gamma_values)
    for j = 1:length(beta_values)
        beta = Beta(i, j);
        gamma = Gamma(i, j);
        
        % Compute R0
        R0 = beta / gamma;
        
        % Compute I* for endemic equilibrium
        if R0 > 1
            I_star(i, j) = N * (1 - 1/R0) * kappa;
        else
            I_star(i, j) = 0;  % Disease-free equilibrium (I* = 0)
        end
    end
end

% Plot the bifurcation diagram
figure;
contourf(Beta, Gamma, I_star, 20, 'LineStyle', 'none'); % Contour plot
colorbar;
xlabel('\beta (Transmission Rate)');
ylabel('\gamma (Recovery Rate)');
title('Bifurcation Diagram of SEIR Model');
set(gca, 'FontSize', 12);
grid on;
hold on;
plot(beta_values, beta_values, 'r--', 'LineWidth', 2);
% annotation('textbox', [0.5, 0.6, 0.075, 0.05], 'String', '$\beta = \gamma$', 'Interpreter', 'latex', ...
%            'FontSize', 12, 'EdgeColor', 'none', 'BackgroundColor', 'white', 'Color', 'red');



% clc; clear; close all;

% Parameter range for beta (bifurcation parameter)
beta_values = linspace(0, 3, 100);  % Varying beta from 0 to 3

% Fixed parameters
N = 1;      % Normalized total population
gamma = 1;  % Recovery rate (fix to 1 for simplicity)
kappa = 1;  % Incubation rate (fix for demonstration)

% Compute R0 for each beta
R0_values = beta_values / gamma;

% Compute endemic equilibrium values of I*
I_star = N * (1 - 1 ./ R0_values) * kappa;  

% Ensure I* is 0 when R0 < 1 (no endemic equilibrium)
I_star(R0_values < 1) = 0;

% Plot bifurcation diagram
figure;
plot(beta_values, I_star, 'b-', 'LineWidth', 2); hold on;
xline(1,'r--', 'LineWidth',1.5);
plot(gamma, 0, 'ro', 'MarkerFaceColor', 'r'); % Bifurcation point

xlabel('R_0');
ylabel('I^* (Steady-State Infected Population)');
title('Transcritical Bifurcation in SEIR Model');
grid on;
set(gca, 'FontSize', 12);

% Annotate the bifurcation point
text(gamma + 0.05, 0.05, 'Bifurcation at \beta = \gamma', 'Color', 'r');


