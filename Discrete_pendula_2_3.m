%Введем матрицы линейной непрерывной системы
A_t = [0, 0, 1, 0;
0, 0, 0, 1;
74.4982, -38.2565, -56.0460, 86.3428;
-114.7697, 101.3459, 86.3428, -228.7317];

B_t = [0; 0; 112.0919; -172.6857];
C = [1, 1, 1, 1];

%Введем параметр шага дискретной системы, построим временное разбиение

h = 0.001;
t_span = (0 : h : 0.09);

% Получаем матрицы дискретной системы
B = integral(@(t)(expm(t * A_t) * B_t),0, h,'ArrayValued',true); 

A = expm(A_t * h);
eig_s = eig(A)
%Получим управление для дискретной системы с помощью формулы Аккермана
%0.1i, -0.1i
%eig_s(1), eig_s(2)
U = -place(A, B, [0.7672, 0.9788, 0.1i,-0.1i])
%U = -place(A, B, [0.7672, 0.9788, 0.98, 0.99])
K = place(A', C', [0.7672, 0.9788, 0.1i,-0.1i]);
%U = -place(A, B, [exp(-265.06819578), exp(-21.41386813), exp(-1), exp(-2)])
%K = place(A', C', [-265.06819578, -21.41386813, -1, -2]);
L = K'


ksi_data = [];
start = [0.1; 0.1; 0; 0];
x_k = start;
ksi_k = [0.; 0.; 0; 0];
% for i = t_span
%     x_k1 = A * x_k + B * U * ksi_k;
%     ksi_k1 = A*ksi_k + B*U*ksi_k + L*C*(x_k - ksi_k);
%     x_data = [x_data, x_k1];
%     x_k = x_k1;
%     ksi_data = [ksi_data, ksi_k1];
%     ksi_k = ksi_k1;
%     
% end
x_data =[];
start = [0.1; 0.1; 0; 0];
x_k = start;
for i = t_span
    x_k1 = A * x_k + B * U * x_k;
    x_data = [x_data, x_k1];
    x_k = x_k1;
end
%Построим графики 

%[ttiks, X] = ode45(@(t,x)(A*x + B*U*x), (0 : h : 10), [0.1, 0.1, 0., 0.]);
%count = length(ksi_data);
count = length(x_data);
utiks = zeros(count, 1);
for k = 1 : count
    utiks(k) = U * (x_data(:, k));
end
%err = C*(ksi_data - x_data);

fileID = fopen('control.txt','w');
fprintf(fileID,'%f\n',utiks);
fclose(fileID);



ttiks = t_span;


subplot(5, 1, 1)
%, ttiks, ksi_data(1,:), 'r'
    plot(ttiks, x_data(1,:), 'b')
    xlabel('t')
    ylabel('\theta_1')
    grid on
    title({'���������� ������� � ����������� � ���� �������� �������� ����� �� ���������.'})
subplot(5, 1, 2)
    plot(ttiks, x_data(2,:), 'b')
    xlabel('t')
    ylabel('\theta_2')
subplot(5, 1, 3)
    plot(ttiks, x_data(3,:), 'b')
    xlabel('t')
    ylabel({'$\dot{\theta}_1$'},'Interpreter','latex')
subplot(5, 1, 4)
    plot(ttiks, x_data(4,:), 'b')
    xlabel('t')
    ylabel({'$\dot{\theta}_2$'},'Interpreter','latex')
subplot(5, 1, 5)
    plot(ttiks, utiks, 'b')
    xlabel('t')
    ylabel('u') 



%{
%Линеаризованная система с кусочно-постоянным управлением из дискретной
%системы
time = 0;
x_0 = start;
t_data = [];
x_data = [];
control = [];
for i = 1 : length(utiks)
    tspan = [time, time + h];
    [t,x] = ode45(@(t,x)(A_t * x + B_t * utiks(i)), tspan, x_0);
    t = t(1: end - 1);
    t_data = [t_data, t'];
    x_0 = x(end, :);
    x_to = x;
    x_to(end, :) = [];
    x_data = [x_data; x_to];
    control = [control, utiks(i) * ones(1, length(x_to))];
    time = time + h;
end
t_data(end)
length(x_data)
length(t_data)
length(control)
subplot(5, 1, 1)
    plot(t_data, x_data(:,1))
    xlabel('t')
    ylabel('\theta_1')
    grid on
    title({'Линеаризованная система'; ' с кусочно-постоянным управлением'})
subplot(5, 1, 2)
    plot(t_data, x_data(:,2))
    xlabel('t')
    ylabel('\theta_2')
subplot(5, 1, 3)
    plot(t_data, x_data(:,3))
    xlabel('t')
    ylabel({'$\dot{\theta}_1$'},'Interpreter','latex')
subplot(5, 1, 4)
    plot(t_data, x_data(:,4))
    xlabel('t')
    ylabel({'$\dot{\theta}_2$'},'Interpreter','latex')
subplot(5, 1, 5)
    plot(t_data, control)
    xlabel('t')
    ylabel('u')
%}