t_k = linspace(0, 10, 1001);
h = 0.01;
A = exp(h);
B = exp(h)-1;
x = [];
x_k = 1; 
u = [];
u_k = 0;

for k = t_k
    theta = -3;
    x_k = A*x_k + B*u_k;
    x = [x, x_k];
    u_k = theta*x_k;
    u = [u, u_k];
end
subplot(2,1,1)
    plot(t_k, x)
    title('График зависимости угла отклонения судна от времени при h = 0.01, \theta = -3')
    xlabel('t_k')
    ylabel('x_k')
    grid on



h = 0.01;
A =[0, exp(h); 0, exp(h)];
B = [(-exp(h)+1);(-exp(h)+1)];
x = [];
x_k = [1;1]; 
u = [];
u_k = 0;
for k = t_k
    theta = -0.05;
    x_k = A*x_k + B*u_k;
    x = [x, x_k];
    u_k = theta*x_k;
    u = [u, u_k];
end
subplot(2,1,2)
    plot(t_k, x)
    title('График зависимости угла отклонения судна от времени при h = 0.01, \theta = -0.05')
    xlabel('t_k')
    ylabel('x_k')
    grid on

