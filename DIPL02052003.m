T = 10;
syms alpha
y = (-8/alpha^5 + 2 * T/alpha^3 - 4 * T^2/ alpha^3 + T^3/alpha^2 + 6 * T/alpha^4) * exp(alpha * T) + 8 / alpha^5;
fy = sqrt(2 * (exp(alpha * T) - 1)/alpha^4 + T * (T / alpha^2 - 2/ alpha^3) * exp(alpha * T));  
res = vpasolve(y,alpha);
f = matlabFunction(fy);
y = f(res)
% x = linspace(0,1);
% y_2 = f(x);
% plot(x,y_2)
% xlabel('$\alpha$','Interpreter','latex')
% ylabel('$f(\alpha)$','Interpreter','latex')