function HomeProject()
A = [0, 1, 0; 1, 0, 1; 0, -1, -7.5];
B = [0; 0; 1];
C = [0, 0, 1];
K = place(A', C', [-10, -9, -11]);
L = K';
H = -place(A, B, [-1, -2, -3]);
G = [A, B*H; L*C, A + B * H - L * C];
a = 7.5;
[ttiks, xticks] = ode45(@(t,x)(G*x), [0, 20], [0.5; 0; 0; 0; 0; 0]);
count = length(ttiks);
zticks = zeros(count, 1);
uticks = zeros(count, 1);
for k = 1 : count
    zticks(k) = C * (xticks(k, 1:3) - xticks(k, 4:6))';
    uticks(k) = H * (xticks(k, 4:6)');
end
% x1 = [];
% x2 = [];
% x3 = [];

% T = linspace(0,20,count);
% for j = 1 : count
%     
%     x1.append(cur(1));
%     x2.append(cur(2));
%     x3.append(cur(3));
% end
subplot(5, 1, 1)
    plot(ttiks, xticks(:, 1), 'b', ttiks,xticks(:, 4), 'r')
subplot(5, 1, 2)
    plot(ttiks, xticks(:, 2), 'b', ttiks,xticks(:, 5), 'r')
subplot(5, 1, 3)
    plot(ttiks, xticks(:, 3), 'b', ttiks,xticks(:, 6), 'r')
subplot(5, 1, 4)
    plot(ttiks, uticks)
[tttiks, xxticks] = ode45(@(t,X)matrr(t,X,a,H),[0,5],[0.5;0;0]);
subplot(5, 1, 5)
    plot(tttiks, xxticks,'b')   
eig(A - L * C)
eig(A + B * H)
end


function y = matrr(t,X,a,H)
  u = H*X;  
  y = [X(2), 0.5*(((1 + X(3))^2)/((1 - X(1))^2) - 1),-(1 + X(3))*X(2)/(1 - X(1)) - a*(1-X(1))*X(3) + (1 - X(1))*u];
end

