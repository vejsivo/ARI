%% zapis matic
A = [0   0   0   1    0    0;
     0   0   0   0    1    0;
     0   0   0   0    0    1;
     7.3809  0   0   0    2    0;
     0  -2.1904  0  -2    0    0;
     0   0  -3.1904  0    0    0];

b1 = [0; 0; 0; 1; 0; 0];
b2 = [0; 0; 0; 0; 1; 0];
b3 = [0; 0; 0; 0; 0; 1];

%% získání eigen hodnot matice přenosu
eig(A)

%% zjištění hodnosti matic kontrolabitlity pro jednotlivé vstupy
Co1 = ctrb(A, b1);
Co2 = ctrb(A, b2);
Co3 = ctrb(A, b3);

n = size(A,1);

disp([rank(Co1) == n, rank(Co2) == n, rank(Co3) == n])