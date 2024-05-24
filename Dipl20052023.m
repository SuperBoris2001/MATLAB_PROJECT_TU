A = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
g = 9.8;
B = [0; -2*g; 0; -g];
C_A_B = [B A*B A*A*B A*A*A*B];
rank(C_A_B)