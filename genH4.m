delta = 0.001;
count = 100;
I = eye(2);
A = [0, 1;
     0, 0];
h = 0.01;
% betta = 0.1;
% A = [0, 1;
%     -1, -2*betta];
Bu = [0;
      1];
Bv = [0;
      1];
C1 = [1, 0];

nx = size(A, 1);
nu = size(Bu, 2);
nw = size(Bv, 2);
nz = size(C1, 1);



alpha = 0.5; % задали  альфа   

lambda = 1; % задали лямбда
      
warning('off', 'CVX:StrictInequalities');
cvx_begin sdp
cvx_solver sdpt3
cvx_precision best
    
variable Y(nx, nx, count) symmetric;
variable Z(nu, nx, count);
variable gammasquare;

minimize(gammasquare);

    for j = 1 : count - 1
        Y(:, :, j + 1) - Y(:, :, j) - h * (A * Y(:, :, j) + Y(:, :, j) * A' + Bv * Bv') > 0
            [Y(:, :, j),       Y(:, :, j) * C1'; 
             C1 * Y(:, :, j),  gammasquare * eye(nz)] > 0        
            Y(:, :, j) > 0;
        [Y(:, :, j),  Z(:, :, j)';
         Z(:, :, j),  (1 - alpha)^2 * gammasquare * eye(nz)] > 0
        Y(:, :, j) > 0;
    end
    [Y(:, :, count),       Y(:, :, count) * C1'; 
         C1 * Y(:, :, count),  gammasquare * eye(nz)] > 0
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
       Ac(:, :, k) = A + Bv * G(:, :, k); 
end 

a = 0.01;

b = 0.99;

c = (a + b)/2;
while(abs(b - a)>delta)
    x = (a + c)/ 2;
    y = (c + b)/ 2;
    
    Gx = max(lambda*funcJ(Ac,Bv, C1, h,count, x, I),(1 - lambda) * funcJ2(Ac, Bv, C1, h, G, x, I))/sqrt(x);% alpha = 0.5
    Gc = max(lambda*funcJ(Ac,Bv, C1, h,count, c, I),(1 - lambda) * funcJ2(Ac, Bv, C1, h, G, c, I))/sqrt(c);% alpha = 0.5
    Gy = max(lambda*funcJ(Ac,Bv, C1, h,count, y, I),(1 - lambda) * funcJ2(Ac, Bv, C1, h, G, y, I))/sqrt(y);% alpha = 0.5
    if (Gc <= Gx)&&(Gc <= Gy)
        a = x;
        b = y;
    end
    if (Gx <= Gc)&&(Gx <= Gy)
        b = c;
        c = x;
    end
    if (Gy <= Gx)&&(Gy <= Gc)
        a = c;
        c = y;
    end
end

function res = funcJ(Ac, Bv, C1, h, count, x_c_y, I)

    %
    nx = size(Ac, 1);
    nz = size(C1, 1);
%    count = size(Ac, 3);

    warning('off', 'CVX:StrictInequalities');
    cvx_begin sdp
        cvx_solver sdpt3
        cvx_precision best
    
        variable Y(nx, nx, count) symmetric; % ошибка
        variable gammasquare;

        minimize(gammasquare);

        Y(:, :, 1) == 0
        for i = 1 : count - 1
            Y(:, :, i + 1) - Y(:, :, i) - h * ((Ac(:, :, i) + 0.5 * x_c_y * I) * Y(:, :, i) + Y(:, :, i) * (Ac(:, :, i) + 0.5 * x_c_y * I)' + Bv * Bv') > 0
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

function res = funcJ2(Ac, Bv, C1, h, G, x_c_y, I)

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
            Y(:, :, i + 1) - Y(:, :, i) + h * ((Ac(:, :, i) + 0.5 * x_c_y * I) * Y(:, :, i) + Y(:, :, i) * (Ac(:, :, i) + 0.5 * x_c_y * I)' + Bv * Bv') < 0
            [Y(:, :, i),  (G(:, :, i)*Y(:, :, i))';
            G(:, :, i)*Y(:, :, i),  gammasquare * eye(nz)] > 0
            Y(:, :, i) > 0;
        end
        [Y(:, :, count), (G(:, :, count)*Y(:, :, count))';
        G(:, :, count)*Y(:, :, count),  gammasquare * eye(nz)] > 0
        Y(:, :, count) > 0
    cvx_end
    warning('on', 'CVX:StrictInequalities');

    res = sqrt(gammasquare);

end
