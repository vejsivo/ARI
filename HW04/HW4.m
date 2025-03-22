s = tf('s');
G = 1 / ((s + 8) * (s + 54));
rltool(G);

%%

[num, den] = tfdata(C_exp, 'v');
Kp = num(1);
Ki = num(2);

disp('Proportional Gain (Kp):'), disp(Kp)
disp('Integral Gain (Ki):'), disp(Ki)