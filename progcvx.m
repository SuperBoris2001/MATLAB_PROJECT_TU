% 1 часть

% A = [-1, 2; 0, -2];
% epsilon = 0.0001;
% I = eye(2);
% cvx_begin sdp
%     variable X(2, 2) symmetric
%     A' * X + X * A <= 0
%     epsilon * I <= X
% cvx_end

% 2 часть
% epsilon = 0.01;
% E = eye(3);
% A = [-4.2386, -0.2026, 0.7193; 2.6649, -2.8342, 0.0175; 0.0344, 0.0005, -3.1772];
% L = [3, 3, 0; 3, 3, 0; 0, 0, -4];
% M = [0, 0, 0; 1, 0, 0; 0, 0, -1];
% I = [0, 0, 0; 1, 0, 0; 0, 0, -1];
% cvx_begin sdp
%     variable X(3, 3) symmetric
%     A' * X + X * A <= 0
%     kron(L, X) + kron(M, (A* X)) + kron(M', (A*X)') <= 0
%     X >= epsilon * E
% cvx_end
% X

% 3 часть
epsilon = 0.01;
E = eye(3);
A = [-1, 2, 0; 0, 0.5, 0; 1, 0, -2];
B = [1, 0; 0, 1; 1, 0];
cvx_begin sdp
    variable Y(3, 3) symmetric
    variable Z(2, 3)
    A*Y + Y*A' + Z'*B' + B*Z <= 0
    Y >= epsilon * E
cvx_end
Y
