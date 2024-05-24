function BallControl1
    g = 9.8;
    mu = 1;
    ms = 1;
    l = 1;
    r = 1;
    s = 1;
    m = 1;
    % Z2 = r * s^2/(m * r * s^2 + 4 * m * s);

%     A = [ 0, 0, 0, 0, 1, 0, 0, 0;
%           0, 0, 0, 0, 0, 1, 0, 0;
%           0, 0, 0, 0, 0, 0, 1, 0;
%           0, 0, 0, 0, 0, 0, 0, 1;
%           0, 15 * g * mu / (l * (28 * ms + 5 * mu)), 0, 0, 0, 0, 0, 0;
%           0, 0, 15 * g * mu / (l * (28 * ms + 5 * mu)), 0, 0, 0, 0, 0;
%           0, - g * mu * (21 * ms + 15 * mu) / (l * (28 * ms * mu + 5 * mu ^ 2)), 0, 0, 0, 0, 0, 0;
%           0, 0, - g * mu * (21 * ms + 15 * mu) / (l * (28 * ms * mu + 5 * mu ^ 2)), 0, 0, 0, 0, 0];
    
    disp(A)
    eig(A)
%     B = [0 0;
%          0 0; 
%          0 0;
%          0 0;
%          4*Z2 0;
%          0 4*Z2;
%          -3*Z2 0;
%          0 -3*Z2];
    B = [0 0;
         0 0;
         0 0;
         0 0;
         20/(28 * ms + 5 * mu) 0;
         0 20/(28 * ms + 5 * mu);
         -15/(28 * ms + 5 * mu) 0;
         0 -15/(28 * ms + 5 * mu)];
    disp(B)
    rank(ctrb(A,B))
    nx = size(A, 1);
    nu = size(B, 2);

    % отключаем предупреждения
    warning('off', 'CVX:StrictInequalities');
    % решаем матричные неравенства    
    cvx_begin sdp
        cvx_precision best
        %cvx_quiet true
        cvx_solver sdpt3
    
        % задаем переменные в неравенствах
        variable Y(nx, nx) symmetric
        variable Z(nu, nx)
            
        % описание ограничений
        Y * A' + A * Y + Z' * B' + B * Z < -0.01 * eye(nx)
        Y > 0
    cvx_end
    % включаем предупреждения
    warning('on', 'CVX:StrictInequalities');

    % выводим найденную матрицу и регулятор
    Y
    G = Z * inv(Y)
    
    
    [ttiks, xticks] = ode45(@(t,x)((A+B*G)*x), (0: 0.01: 5), [0, 0, 0.1, 0.1, 0, 0, 0, 0]);
    
    subplot(8, 1, 1)
    plot(ttiks,xticks(:, 1),'b')
    xlabel('t')
    ylabel('x')
    subplot(8, 1, 2)
    plot(ttiks,xticks(:, 2),'b')
    xlabel('t')
    ylabel('y')
    subplot(8, 1, 3)
    plot(ttiks,xticks(:, 3),'b')
    xlabel('t')
    ylabel('X')
    subplot(8, 1, 4)
    plot(ttiks,xticks(:, 4),'b')
    xlabel('t')
    ylabel('Y')
 
    
% 4*mu*Derivative(X(t), (t, 2))/3 + mu*Derivative(x(t), (t, 2)) + (-g*mu/sqrt(l**2 - X(t)**2 - Y(t)**2) - 4*(X(t)*Derivative(X(t), t) + Y(t)*Derivative(Y(t), t))**2/(3*(l**2 - X(t)**2 - Y(t)**2)**2) + 4*(X(t)*Derivative(X(t), t) + Y(t)*Derivative(Y(t), t))*(2*X(t)*Derivative(X(t), t) + 2*Y(t)*Derivative(Y(t), t))/(3*(l**2 - X(t)**2 - Y(t)**2)**2) + 4*(X(t)*Derivative(X(t), (t, 2)) + Y(t)*Derivative(Y(t), (t, 2)) + Derivative(X(t), t)**2 + Derivative(Y(t), t)**2)/(3*(l**2 - X(t)**2 - Y(t)**2)))*X(t)    

% 4*mu*Derivative(Y(t), (t, 2))/3 + mu*Derivative(y(t), (t, 2)) + (-g*mu/sqrt(l**2 - X(t)**2 - Y(t)**2) - 4*(X(t)*Derivative(X(t), t) + Y(t)*Derivative(Y(t), t))**2/(3*(l**2 - X(t)**2 - Y(t)**2)**2) + 4*(X(t)*Derivative(X(t), t) + Y(t)*Derivative(Y(t), t))*(2*X(t)*Derivative(X(t), t) + 2*Y(t)*Derivative(Y(t), t))/(3*(l**2 - X(t)**2 - Y(t)**2)**2) + 4*(X(t)*Derivative(X(t), (t, 2)) + Y(t)*Derivative(Y(t), (t, 2)) + Derivative(X(t), t)**2 + Derivative(Y(t), t)**2)/(3*(l**2 - X(t)**2 - Y(t)**2)))*Y(t)

%     subplot(8, 1, 5)
%     plot(ttiks,xticks(:, 5),'b')
%     
%     plot(ttiks,xticks(:, 6),'b');
%     plot(ttiks,xticks(:, 7),'b');
%     plot(ttiks,xticks(:, 8),'b');
    % для проверки находим собственные числа
    eig(A + B * G)
%     n = size(A, 1)
%     % nx = size(B, 1)
%     warning('off', 'CVX:StrictInequalities');
%     cvx_begin sdp
%         cvx_solver sdpt3
%         cvx_precision best
%         cvx_quiet true
%         variable Y(n, n) symmetric
%         
%     cvx_end
end

function res = matrixA2(x, g, mu, l)
    res = [ 0, 0, 0, 0, 1, 0, 0, 0;
          0, 0, 0, 0, 0, 1, 0, 0;
          0, 0, 0, 0, 0, 0, 1, 0;
          0, 0, 0, 0, 0, 0, 0, 1;
          0, 15 * g * mu / (l * (28 * ms + 5 * mu)), 0, 0, 0, 0, 0, 0;
          0, 0, 15 * g * mu / (l * (28 * ms + 5 * mu)), 0, 0, 0, 0, 0;
          0, gamma1(x[2], x[3], x[6], x[7], g, mu, l), 0, 0, 0, 0, 0, 0;
          0, 0, - g * mu * (21 * ms + 15 * mu) / (l * (28 * ms * mu + 5 * mu ^ 2)), 0, 0, 0, 0, 0];
end

function res = gamma1(X, Y, dX, dY, g, mu, l)
    res = -g*mu/sqrt(l^2 - X^2 - Y^2)
    - 4*(X*dX + Y*dY^2)/(3*(l^2 - X^2 - Y^2)^2)
    + 4*(X*dX + Y*dY*(2*X*dX + 2*Y*dY))/(3*(l^2 - X^2 - Y^2)^2)
    + 4*(X*dX + Y*dY + dX^2 + dY^2)/(3*(l^2 - X^2 - Y^2)) 
end