% A = [0, 1;
%      0, 0];
% betta = 0.1;
% A = [0, 1;
%     -1, -2*betta];
k = 1;
m = 1;
A = [0 1 0 0;
     -2*k/m 0 k/m 0;
     0 0 0 1;
     k/m 0 -k/m 0];
Bu = [0;
      1;
      0;
      1];
Bv = [0;
      1;
      0;
      1];
C1 = [-1, 0, 1, 0];


% определяем размеры матриц
nx = size(A, 1);
nu = size(Bu, 2);
nw = size(Bv, 2);
nz = size(C1, 1);

% кол-во узлов в сетке
count = 100;

alpha = 0.1;

tticks = (0:0.01:0.99);

% формируем сетку alpha
alpha_count = 0;

alpha_min = 0.02;

alpha_max = 0.98;

h_a = 0.01; % задаём шаг \alpha
alpha_array = [];

alpha_flag = alpha_min;

% формируем файл из коэффициентов линейной обратной связи 
filename = 'genH2mat2_theta_3.dat';
        fid = fopen(filename, 'wt');
        fprintf(fid, 'alpha,J1,J2,theta_1,theta_2,theta_3,theta_4,tticks\n');  % header
        fclose(fid);

while alpha_flag <= alpha_max
    alpha_count = alpha_count + 1;
    if alpha_flag < 0.05
        alpha_array = [alpha_array, alpha_flag];
        h_a = 0.01;
        alpha_flag = alpha_flag + h_a;
    elseif alpha_flag <= 0.95
        alpha_array = [alpha_array, alpha_flag];
        h_a = 0.05;
        alpha_flag = alpha_flag + h_a;
    elseif alpha_flag <= 0.99
        alpha_array = [alpha_array, alpha_flag];
        h_a = 0.01;
        alpha_flag = alpha_flag + h_a;
    end
end
alpha_count = alpha_count + 1;
alpha_array = [alpha_array, alpha_flag]; % сетка параметра alpha

length_alpha_array = size(alpha_array, 2);

% находим зависимость J1 от J2 при изменении alpha 

theta_array = zeros(2,length_alpha_array); % список J

for i = 1:length_alpha_array
    
    alpha = alpha_array(i);
    h = tticks(1,2) - tticks(1,1);
    
    theta_1 = zeros(1,count); % список Theta1
    theta_2 = zeros(1,count); % список Theta2
    theta_3 = zeros(1,count); % список Theta3
    theta_4 = zeros(1,count); % список Theta4
    warning('off', 'CVX:StrictInequalities');
    cvx_begin sdp
      cvx_solver sdpt3
      cvx_precision best
    
    variable Y(nx, nx, count) symmetric;
    variable Z(nu, nx, count);
    variable gammasquare;

    minimize(gammasquare);

    Y(:, :, 1) == eye(nx, nx);
    for j = 1 : count - 1
        Y(:, :, j + 1) - Y(:, :, j) - h * (A * Y(:, :, j) + Bu * Z(:, :, j) + Y(:, :, j) * A' + Z(:, :, j)' * Bu' + Bv * Bv') == 0
        [Y(:, :, j),       Y(:, :, j)' * C1'; 
         C1 * Y(:, :, j),  alpha^2 * gammasquare * eye(nz)] > 0
        [Y(:, :, j),  Z(:, :, j)';
         Z(:, :, j),  (1 - alpha)^2 * gammasquare * eye(nz)] > 0
        Y(:, :, j) > 0;
    end
    [Y(:, :, count),       Y(:, :, count)' * C1'; 
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
       theta_1(1,k) = G(1,1,k);
       theta_2(1,k) = G(1,2,k);
       theta_3(1,k) = G(1,3,k);
       theta_4(1,k) = G(1,4,k);
    end 
    
    J1 = funcJ1(Ac, Bv, C1, h);
    J2 = funcJ2(Ac, Bv, C1, h, G);
    alpha_e = alpha*ones(1,count);
    J1_e = J1*ones(1,count);
    J2_e = J2*ones(1,count);
    
    theta_array(1,i) = J1;
    theta_array(2,i) = J2;
    dlmwrite(filename, [alpha_e;J1_e;J2_e;theta_1;theta_2;theta_3;theta_4;tticks]', ...
            'delimiter', ',', 'precision', '%.6e', '-append');
end

ax1 = subplot(1,1,1);
plot(theta_array(1,1:end), theta_array(2,1:end));
xlabel(ax1, 'J_1');
ylabel(ax1, 'J_2');


% функции вычисления J_1 и J_2
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
            [Y(:, :, i),       Y(:, :, i)' * C1'; 
             C1 * Y(:, :, i),  gammasquare * eye(nz)] > 0
            Y(:, :, i) > 0;
        end
        [Y(:, :, count),       Y(:, :, count)' * C1'; 
         C1 * Y(:, :, count),  gammasquare * eye(nz)] > 0
        Y(:, :, count) > 0
    cvx_end
    warning('on', 'CVX:StrictInequalities');

    res = sqrt(gammasquare);

end



function res = funcJ2(Ac, Bv, C1, h, G)

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