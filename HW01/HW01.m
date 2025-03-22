%% ploting the simulink model
% Extract time and values

time = simout1.Time;
values1 = simout1.Data;
values2 = simout2.Data;

figure(1)
plot(time, values1, 'b', 'LineWidth', 0.6) % First signal (blue)
hold on
plot(time, values2, 'r', 'LineWidth', 0.6) % Second signal (red)
xlabel('Time (s)')
ylabel('Amplitude');
title('Time response');
grid on;
legend('original system', 'linearized system'); % Adjust if needed
print(1, 'out1', '-depsc');