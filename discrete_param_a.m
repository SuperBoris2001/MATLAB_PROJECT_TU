%Линеаризованная система
a = 7.5;
A = [0, 1, 0; 1, 0, 1; 0, -1, -a];
B = [0; 0; 1];
C = [0, 0, 1];

K = place(A', C', [-1, -2, -3]);
L = K';

H = -place(A, B, [-3, -2+1i, -2-1i]);

%Полная система
[ttiks, X] = ode45(@(t, X) (to_ode_1(t, X, a, H)), (0: 0.001: 1), [0.5, 0, 0]);
uticks = [];
%Управление в виде функции Хевисайда
for k = (0: 0.001: 1)
    uticks = [uticks; heaviside(k-0.5)];
end
%Стабилизирующее управление theta*X
%{
for k =(1:length(ttiks))
   uticks = [uticks; H*X(k,:)'];
end
%}

subplot(4, 1, 1)
    plot(ttiks, X(:, 1), 'b')
    xlabel('t')
    ylabel('x')
    grid on;
title({'Дискретная нелинейная система'})
subplot(4, 1, 2)
    plot(ttiks, X(:, 2), 'b')
    xlabel('t')
    ylabel('x`')
    grid on;
subplot(4, 1, 3)
    plot(ttiks, X(:, 3), 'b')
    xlabel('t')
    ylabel('I')
    grid on;
subplot(4, 1, 4)
    plot(ttiks, uticks, 'r')
    xlabel('t')
    ylabel('u')
    grid on;

a_data = zeros(length(ttiks)-1,1);
y_data = [[],[],[]];
C_data = [[],[],[]];
H_data = zeros(length(ttiks)-1,1);
T_data = zeros(length(ttiks)-1,1);
h=0.001;
for k = (1:length(ttiks)-1)
    y_data(k,1) = X(k+1,1) - X(k,1) - X(k,2)*h;
    y_data(k,2) = X(k+1,2) - X(k,2) - 0.5*((1+X(k,3)^2)/(1-X(k,1))^2 - 1)*h;
    y_data(k,3) = X(k+1,3) - X(k,3) + ((1+X(k,3))/(1-X(k,1))*X(k,2) - (1 - X(k,1))*uticks(k))*h;
    C_data(k,1) = 0;
    C_data(k,2) = 0;
    C_data(k,3) = -h*X(k,3)*(1 - X(k,1));
    for j = (1:k)
        H_data(k) = H_data(k) + C_data(j,:)*C_data(j,:)';
    end
end
for k = (1:length(ttiks)-2) 
    a_data(k+1) = a_data(k) + inv(H_data(k+1) + C_data(k+1,:)*C_data(k+1,:)')*C_data(k+1,:)*(y_data(k+1,:) - C_data(k+1,:)*a_data(k))';
end
a_data(end)

function dxdt = to_ode_1(t, X, a, H)
    u = heaviside(t-0.5);
    %u = H*X;
    x_1 = X(1);
    x_2 = X(2);
    x_3 = X(3);
    dx_1dt = x_2;
    dx_2dt = 1/2*(((1+x_3)^2)/((1-x_1)^2)-1);
    dx_3dt = -x_2*(1+x_3)/(1-x_1) - a*(1-x_1)*x_3 + (1-x_1)*u;
    dxdt = [dx_1dt; dx_2dt; dx_3dt];
end
