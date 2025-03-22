% Define Kp and Ki values for the stable region
Kp = [-2, 0, 10, 10, -2]; 
Ki = [0, 0, 10, 10, 0]; 

% Create figure
figure;
hold on;
grid on;

% Plot the stable region
fill(Kp, Ki, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Labels and title
xlabel('$K_p$', 'Interpreter', 'latex');
ylabel('$K_I$', 'Interpreter', 'latex');
title('Oblast hodnot zesílení se stabilní uzavřenou smyčkou', 'Interpreter', 'latex');

% Legend
legend({'Stabilní oblast'}, 'Interpreter', 'latex', 'Location', 'northeast');

% Axis limits
xlim([-2, 10]);
ylim([0, 10]);

% Add LaTeX caption using annotation
annotation('textbox', [0.35, -0.12, 0.3, 0.05], 'String', '\textbf{Obrázek 1:} Oblast zesílení kaskádního kompenzátoru', ...
    'Interpreter', 'latex', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

hold off;
