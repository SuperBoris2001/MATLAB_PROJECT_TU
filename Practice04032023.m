%A = [0, 1, 0; 1, 0, 1; 0, -1, -a];
%B = [0; 0; 1];
A = [1, 1; 0, 2];
B = [0; 1];
C = [1, 0];
L = [4; 7.25];
H = [-5, -5];
%G = [A, zeros(2);L*C,A-L*C]
G = [A, B*H; L*C, A + B * H - L * C];
[ttiks, xticks] = ode45(@(t,x)(G*x), [0, 20], [1; 1; 0; 0]);
count = length(ttiks);
zticks = zeros(count, 1);
for k = 1 : count
 zticks(k) =C * (xticks(k,1:2) - xticks(k,3:4))';
end
subplot(3, 1, 1);
 plot(ttiks,xticks(:, 1),'b', ttiks,xticks(:, 3),'r');
subplot(3, 1, 2);
 plot(ttiks,xticks(:, 2),'b', ttiks,xticks(:, 4),'r');
subplot(3, 1, 3);
 plot(ttiks, zticks);
% ######## ###
eig(A-L*C);
K = place(A',C',[-0.5 + 1i, -0.5 - 1i]);
L = K';
eig(A-L*C);
H = -place(A, B, [-1 + 1i, -1 - 1i]);
eig(A+ B*H);
