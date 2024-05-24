A = [0, 1;
     0, 0];
Bu = [0;
      1];
Bv = [0;
      1];
C1 = [1, 0];

% определяем размеры матриц
nx = size(A, 1);
nu = size(Bu, 2);
nw = size(Bv, 2);
nz = size(C, 1);

% кол-во узлов в сетке
count = 100;

alpha = 0.5;

h = 0.01;

warning('off', 'CVX:StrictInequalities');
cvx_begin sdp
    cvx_solver sdpt3
    cvx_precision best
    
    variable Y(nx, nx, count) symmetric;
    variable Z(nu, nx, count);
    variable gammasquare;

    minimize(gammasquare);

    for i = 1 : count - 1
        Y(:, :, i + 1) - Y(:, :, i) + h * (A * Y(:, :, i) + Bu * Z(:, :, i) + Y(:, :, i) * A' + Z(:, :, i)' * Bu' + Bv * Bv') < 0
        [Y(:, :, i),       Y(:, :, i) * C1'; 
         C1 * Y(:, :, i),  alpha^2 * gammasquare * eye(nz)] > 0
        [Y(:, :, i),  Z(:, :, i)';
         Z(:, :, i),  (1 - alpha)^2 * gammasquare * eye(nz)] > 0
        Y(:, :, i) > 0;
    end
    [Y(:, :, count),       Y(:, :, count) * C1'; 
     C1 * Y(:, :, count),  alpha^2 * gammasquare * eye(nz)] > 0
    [Y(:, :, count),  Z(:, :, count)';
     Z(:, :, count),  (1 - alpha)^2 * gammasquare * eye(nz)] > 0
    Y(:, :, count) > 0
cvx_end
warning('on', 'CVX:StrictInequalities');

% формируем матрицы замкнутой системы
Ac = zeros(nx, nx, count);
G = zeros(nu, nx, count);
for k = 1 : count
   G(:, :, k) = Z(:, :, k) * inv(Y(:, :, k));
   Ac(:, :, k) = A + Bu * G(:, :, k); 
end


res = funcJ1(Ac, Bv, C1, h)


warning('off', 'CVX:StrictInequalities');
cvx_begin sdp
    cvx_solver sdpt3
    cvx_precision best
    
    variable Y(nx, nx, count) symmetric;
    variable gammasquare;

    minimize(gammasquare);

    for i = 1 : count - 1
        Y(:, :, i + 1) - Y(:, :, i) + h * (Ac(:, :, i) * Y(:, :, i) + Y(:, :, i) * Ac(:, :, i)' + Bv * Bv') < 0
        [Y(:, :, i),  G(:, :, i)';
         G(:, :, i),  gammasquare * eye(nz)] > 0
        Y(:, :, i) > 0;
    end
    [Y(:, :, i),  G(:, :, i)';
     G(:, :, i),  gammasquare * eye(nz)] > 0
    Y(:, :, count) > 0
cvx_end
warning('on', 'CVX:StrictInequalities');

J2 = sqrt(gammasquare)


function res = funcJ1(Ac, Bv, C1, h)

    %
    nx = size(Ac, 1);
    nz = size(C1, 1);
    count = size(Ac, 3);

    warning('off', 'CVX:StrictInequalities');
    cvx_begin sdp
        cvx_solver sdpt3
        cvx_precision best
    
        variable Y(nx, nx, count) symmetric;
        variable gammasquare;

        minimize(gammasquare);

        for i = 1 : count - 1
            Y(:, :, i + 1) - Y(:, :, i) + h * (Ac(:, :, i) * Y(:, :, i) + Y(:, :, i) * Ac(:, :, i)' + Bv * Bv') < 0
            [Y(:, :, i),       Y(:, :, i) * C1'; 
             C1 * Y(:, :, i),  gammasquare * eye(nz)] > 0
            Y(:, :, i) > 0;
        end
        [Y(:, :, count),       Y(:, :, count) * C1'; 
         C1 * Y(:, :, count),  gammasquare * eye(nz)] > 0
        Y(:, :, count) > 0
    cvx_end
    warning('on', 'CVX:StrictInequalities');

    res = sqrt(gammasquare);

end

