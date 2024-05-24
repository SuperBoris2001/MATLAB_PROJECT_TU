% % A = [0, 1, 0; 0, 0, 1; 10, 1, 5];
% % B = [0; 1; 1];
% % C = [B, A*B, A^2*B];
% 
% A = [2 -1 0 -1; 7 -2 4 -5; -11 4 -6 9; -10 3 -7 9];
% B = [3; 3; -5; -2];
% C = [6 -2 4 -6];
% 
% C_AB = [B, A * B,A^2 * B, A^3 * B]
% rank(C_AB)
% 
% [e,D] = eig(C_AB)
% 
% G_AC = [C; C * A; C * A^2; C * A^3]
% rank(G_AC)
% 
% res = null(G_AC);
% 
% res2 = [3 5 0.4614 0.6222; 3 5 -0.3099 0.5513; -5 -9 0.3099 -0.5513;-2 -4 0.7713 0.0709]
% res3 = null(res2)
% res3(1)
% res3(2)
% res4 = [3; 3; -5; -2]*res3(1) + [5; 5; -9; -4]*res3(2)
% 
% T = [3, 0.3536, 0 ,- 0.4614; 3, 0.3536, 0, 0.3099; -5, -0.3536, 0, -0.3099; -2, 0, 1, -0.7713];
% 
% rank(T)
% Anew = inv(T) * A * T
% % eig(A)
% % 
% % H = -place(A, B, [-1, -0.2554+1.3227i, -0.2554-1.3227i]);
% % 
% % eig(A + B*H)
% 
% dop1 = [1 0 ; -5.186 -1]
% eig(dop1)
% 
% dop2 = [2 0; -2.8281 1]
% eig(dop2)


%roots([1 3 30 2 1])

%roots([1 1 30 10 1])

roots([2 1 10 10 2])

% syms w
% y = (3 - 9 * w^2 + 2 * w^4)/(1 + w^2 + 2 * w^4);
% diff(y)
% y = expand(y)
% y = expand(y)
% 
% x = (5 - 6 * w^2)/(1 + 6 * w^2);
% diff(x)
% x = expand(x)
% 
% yf = y/x
% 
% yf = expand(yf)

A = [0 0 1 0;
    0 0 0 1;
    74.49823574469526 -38.256590310670866 -56.04604080240516 86.34286705234221;
    -114.76977093201259 101.34592586957349 86.34286705234221 -228.73177490708963];
B = [0;
    0;
    112.0919;
    -172.6856];
h = 0.01;

%y = integral(@(tau)(expm(A * tau) * B), 0, h,'ArrayValued',true)

Q = [0 0 0 0;
     0 0 0 0;
     0 0 0 0;
     0 0 0 0]
R = 1;
[X,L,G] = care(A,B,Q,R)

F = [A - B*G]
[ttiks, xticks] = ode45(@(t,x)(F*x), (0: 0.0001: 10), [0.6, 0.1, 0, 0]);

count = length(ttiks);
zticks = zeros(count, 1);
uticks = zeros(count, 1);
for k = 1 : count
uticks(k) = -G * (xticks(k,1:4))';
end


subplot(6, 2, 1)
plot(ttiks,xticks(:, 1),'b')
xlabel('t')
ylabel('\theta_1')
subplot(6, 2, 2)
plot(ttiks,xticks(:, 2),'b')
xlabel('t')
ylabel('\theta_2')
subplot(6, 2, 3)
plot(ttiks, xticks(:, 3),'b')
xlabel('t')
ylabel({'$\dot{\theta}_1$'},'Interpreter','latex')
subplot(6, 2, 4)
plot(ttiks, xticks(:, 4),'b')
xlabel('t')
ylabel({'$\dot{\theta}_2$'},'Interpreter','latex')


% Управление по состоянию наблюдателя
subplot(6,2,5)
plot(ttiks, uticks)
xlabel('t')
ylabel('u')


[eig(A)   eig(A-B*G)]


