A = [0, 1, 0; 1, 0, 1; 0, -1, -7.5];
B = [0; 0; 1];
C = [0, 0, 1];
K = place(A', C', [-10, -9, -11]);
L = K';
H = -place(A, B, [-1, -2, -3]);
G = [A, B*H; L*C, A + B * H - L * C];
[ttiks, xticks] = ode45(@(t,x)(G*x), [0, 20], [0.5; 0; 0; 0; 0; 0]);
count = length(ttiks);
zticks = zeros(count, 1);
uticks = zeros(count, 1);
for k = 1 : count
    zticks(k) = C * (xticks(k, 1:3) - xticks(k, 4:6))';
    uticks(k) = H * (xticks(k, 4:6)');
end
subplot(4, 1, 1)
    plot(ttiks, xticks(:, 1), 'b', ttiks,xticks(:, 4), 'r')
subplot(4, 1, 2)
    plot(ttiks, xticks(:, 2), 'b', ttiks,xticks(:, 5), 'r')
subplot(4, 1, 3)
    plot(ttiks, xticks(:, 3), 'b', ttiks,xticks(:, 6), 'r')
subplot(4, 1, 4)
    plot(ttiks, uticks)
% œ–Œ¬≈– » 

eig(A - L * C)
eig(A + B * H)