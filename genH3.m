delta = 0.001;
count = 100;
I = eye(2);
A = [0, 1;
     0, 0];
h = 0.01;
% betta = 0.1;
% A = [0, 1;
%     -1, -2*betta];
% Bu = [0;
%       1];
Bv = [0;
      1];
C1 = [1, 0];

a = 0.01;

b = 0.99;

c = (a + b)/2;
while(abs(b - a)>delta)
    x = (a + c)/ 2;
    y = (c + b)/ 2;
    Gx = funcJ(A + 0.5 * x * I,Bv, C1, h)/sqrt(x);% alpha = 0.5
    Gc = funcJ(A + 0.5 * c * I,Bv, C1, h)/sqrt(c);% alpha = 0.5
    Gy = funcJ(A + 0.5 * y * I,Bv, C1, h)/sqrt(y);% alpha = 0.5
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

function res = funcJ(Ac, Bv, C1, h)

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

        Y(:, :, 1) == 0
        for i = 1 : count - 1
            Y(:, :, i + 1) - Y(:, :, i) + h * (Ac(:, :, i) * Y(:, :, i) + Y(:, :, i) * Ac(:, :, i)' + Bv * Bv') == 0
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
