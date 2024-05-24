%Линеаризованная система
a = 5;
A = [0, 1, 0; 1, 0, 1; 0, -1, -a];
B = [0; 0; 1];
C = [0, 0, 1];

K = place(A', C', [-1, -2, -3]);
L = K'

H = -place(A, B, [-1, -2+1i, -2-1i])

%Управление по состоянию наблюдателя
G = [A, B*H; L*C, A + B * H - L * C];
[ttiks, X] = ode45(@(t,x)(G*x), [0, 10], [0.25, 0, 0, 0, 0, 0]);
subplot(4, 2, 1)
    plot(ttiks, X(:, 1), 'b', ttiks, X(:, 4), 'r')
    xlabel('t')
    ylabel('x')
    title({'Линеаризованная система'; ' с управлением по состоянию наблюдателя'})
subplot(4, 2, 3)
    plot(ttiks, X(:, 2), 'b', ttiks, X(:, 5), 'r')
    xlabel('t')
    ylabel('x`')
subplot(4, 2, 5)
    plot(ttiks, X(:, 3), 'b', ttiks, X(:, 6), 'r')
    xlabel('t')
    ylabel('I')
count = length(ttiks);
utiks = zeros(count, 1);
for k = 1 : count
    X(k, 4:6);
    utiks(k) = H * (X(k, 4:6)');
end
subplot(4, 1, 4)
    plot(ttiks, utiks)
    xlabel('t')
    ylabel('u')

%Полная система
[ttiks, X_n] = ode45(@(t, X) to_ode_2(t, X, a, A, B, H, L, C), [0, 20], [0.25, 0, 0, 0, 0, 0]);

subplot(4, 2, 2)
    plot(ttiks, X_n(:, 1), 'b', ttiks, X_n(:, 4), 'r')
    xlabel('t')
    ylabel('x')
title({'Нелинейная система '; 'с подставленным управлением'})
subplot(4, 2, 4)
    plot(ttiks, X_n(:, 2), 'b', ttiks, X_n(:, 5), 'r')
    xlabel('t')
    ylabel('x`')
subplot(4, 2, 6)
    plot(ttiks, X_n(:, 3), 'b', ttiks, X_n(:, 6), 'r')
    xlabel('t')
    ylabel('I')


function dxdt = to_ode_1(t,X, a, H)
    u = H*X;
    x_1 = X(1);
    x_2 = X(2);
    x_3 = X(3);
    dx_1dt = x_2;
    dx_2dt = 1/2*(((1+x_3)^2)/((1-x_1)^2)-1);
    dx_3dt = -x_2*(1+x_3)/(1-x_1) - a*(1-x_1)*x_3 + (1-x_1)*u;
    dxdt = [dx_1dt; dx_2dt; dx_3dt];
end

function dxdt = to_ode_2(t, X, a, A, B, H, L, C)
    x_1 = X(1);
    x_2 = X(2);
    x_3 = X(3);
    ksi_1 = X(4);
    ksi_2 = X(5);
    ksi_3 = X(6);
    ksi = [ksi_1; ksi_2; ksi_3];
    dksidt = A*ksi + B*H*ksi + L*C*(X(1:3)-ksi);
    dksi_1dt = dksidt(1);
    dksi_2dt = dksidt(2);
    dksi_3dt = dksidt(3);
    u = H*ksi;
    dx_1dt = x_2;
    dx_2dt = 1/2*(((1+x_3)^2)/((1-x_1)^2)-1);
    dx_3dt = -x_2*(1+x_3)/(1-x_1) - a*(1-x_1)*x_3 + (1-x_1)*u;
    dxdt = [dx_1dt; dx_2dt; dx_3dt; dksi_1dt; dksi_2dt; dksi_3dt];
end
