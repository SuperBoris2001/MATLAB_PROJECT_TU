A = [0 , 1; 0 , 1];
B = [0 ; -1];
ro = 1;
R = ro;
Q = eye(2);
I = eye(2);
ksi0 = [1;1];
rho = 
for k = 0:100
    X = care(A,B,Q,I);
    Ac = A - B*R*B'*X;
    