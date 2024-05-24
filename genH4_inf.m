function genH4

    delta = 0.001;

    lamlist = 0.02:0.02:0.98;
    count = length(lamlist);
    vals_list = zeros(count, 5);
    k = 1;
    filename = 'genH2mat4inf.dat';
        fid = fopen(filename, 'wt');
        fprintf(fid, 'lambda,alpha,G1,J1,G2,J2\n');  % header
        fclose(fid);
    for lam = lamlist
        a = 0.01;
        b = 10.0;
        c = (a + b) / 2;
        while(abs(b - a) > delta)
            x = (a + c) / 2;
            y = (c + b) / 2;
    
            [Gx, Ax, Cx1, Cx2] = bound_funcs(x, lam);
            [Gc, Ac, Cc1, Cc2] = bound_funcs(c, lam);
            [Gy, Ay, Cy1, Cy2] = bound_funcs(y, lam);
            if (Gc <= Gx) && (Gc <= Gy)
                a = x;
                b = y;
            end
            if (Gx <= Gc) && (Gx <= Gy)
                b = c;
                c = x;
            end
            if (Gy <= Gx) && (Gy <= Gc)
                a = c;
                c = y;
            end
        end
        val_list(k, 1) = c;
        % определяем верчнюю границу изменения alpha
        nx = size(Ac, 1);
        alpha = 0.0 : 0.1 : 10;
        tbl = arrayfun(@(x)max(real(eig(Ac + 0.5 * x * eye(nx)))), alpha)
        bmax = max(alpha(tbl < 0));
        % считаем значения 1-го функционала
        a = 0.01;
        b = bmax;
        c = (a + b) / 2;
        while(abs(b - a) > delta)
            x = (a + c) / 2;
            y = (c + b) / 2;
    
            [Gx, Jx] = func_Ja(x, Ac, Cc1);
            [Gc, Jc] = func_Ja(c, Ac, Cc1);
            [Gy, Jy] = func_Ja(y, Ac, Cc1);
            if (Gc <= Gx) && (Gc <= Gy)
                a = x;
                b = y;
            end
            if (Gx <= Gc) && (Gx <= Gy)
                b = c;
                c = x;
            end
            if (Gy <= Gx) && (Gy <= Gc)
                a = c;
                c = y;
            end
        end
        val_list(k, 2) = Gc;
        val_list(k, 3) = Jc;
        
         % считаем значения 2-го функционала
        a = 0.01;
        b = bmax;
        c = (a + b) / 2;
        while(abs(b - a) > delta)
            x = (a + c) / 2;
            y = (c + b) / 2;
    
            [Gx, Jx] = func_Ja(x, Ac, Cc2);
            [Gc, Jc] = func_Ja(c, Ac, Cc2);
            [Gy, Jy] = func_Ja(y, Ac, Cc2);
            if (Gc <= Gx) && (Gc <= Gy)
                a = x;
                b = y;
            end
            if (Gx <= Gc) && (Gx <= Gy)
                b = c;
                c = x;
            end
            if (Gy <= Gx) && (Gy <= Gc)
                a = c;
                c = y;
            end
        end
        val_list(k, 4) = Gc;
        val_list(k, 5) = Jc; 
        fprintf('--------------------------------------------------------------------------------------------------------\n');
        fprintf('lambda is %4.2f | alpha is %4.6f | G1 is %4.6f | J1 is %4.6f | G2 is %4.6f | J2 is %4.6f',lam, val_list(k,1),val_list(k,2),val_list(k,3),val_list(k,4),val_list(k,5))
        fprintf('\n');
              
        dlmwrite(filename, [lam; val_list(k,1);val_list(k,2);val_list(k,3);val_list(k,4);val_list(k,5)]', ...
            'delimiter', ',', 'precision', '%.6e', '-append');
       
        k = k + 1;
        
    end
%     filename = 'genH2mat4inf.dat';
%     fid = fopen(filename, 'wt');
%     fprintf(fid, 'lambda,alpha,G1,J1,G2,J2\n');  % header
%     fclose(fid);      
%     dlmwrite(filename, [lamlist; val_list]', ...
%         'delimiter', ',', 'precision', '%.6e', '-append');
    
    ax1 = subplot(1,1,1);
    plot(val_list(1:end, 2), val_list(1:end, 4));
    xlabel(ax1, '\theta_1');
    ylabel(ax1, '\theta_2');

end


function [val, Ac, Cc1, Cc2] = bound_funcs(alpha, lam)

    % задаем матрицы системы
    A = [0, 1;
         0, 0];
    Bu = [0;
          1];
    Bv = [0;
          1];
    C1 = [1, 0];
    D1 = 0;
    C2 = [0, 0];
    D2 = 1;
    
    % определяем размер
    nx = size(A, 1);
    nu = size(Bu, 2);
    nw = size(Bv, 2);
    nz = size(C1, 1);
       
    %
    [~, G, gval] = genH2(A + 0.5 * alpha * eye(nx), Bu, Bv, C1, D1, C2, D2, lam);
    % формируем матрицы замкнутой системы
    Ac = A + Bu * G; 
    Cc1 = C1 + D1 * G; 
    Cc2 = C2 + D2 * G; 
    val = gval / sqrt(alpha);

end

function [Y, G, gval] = genH2(A, Bu, Bv, C1, D1, C2, D2, lam)

    nx = size(A, 1);
    nu = size(Bu, 2);
    nw = size(Bv, 2);
    nz1 = size(C1, 1);
    nz2 = size(C2, 1);

    warning('off', 'CVX:StrictInequalities');
    cvx_begin sdp
        cvx_solver sdpt3
        cvx_precision best
        cvx_quiet true
    
        variable Y(nx, nx) symmetric
        variable Z(nu, nx)
        variable gammasquare

        minimize(gammasquare)

        A * Y + Y * A' + Bu * Z + Z' * Bu' + Bv * Bv' < 0
        [Y,                         Y * C1' + Z' * D1'; 
         C1 * Y + D1 * Z,  lam^2 * gammasquare * eye(nz1)] > 0        
        [Y,                         Y * C2' + Z' * D2';
         C2 * Y + D2 * Z,  (1 - lam)^2 * gammasquare * eye(nz2)] > 0
        Y > 0
    cvx_end
    warning('on', 'CVX:StrictInequalities');
    
    G = Z * inv(Y);
    gval = sqrt(gammasquare);

end

function [jval, gval] = func_Ja(alpha, Ac, Cc)

    % задаем матрицы системы
    Bu = [0;
          1];
    Bv = [0;
          1];
    
    % определяем размер
    nx = size(Ac, 1);
    nu = size(Bu, 2);
    nw = size(Bv, 2);
    nz = size(Cc, 1);
       
    %
    Aa = Ac + 0.5 * alpha * eye(nx);
    
    warning('off', 'CVX:StrictInequalities');
    cvx_begin sdp
        cvx_solver sdpt3 
        cvx_quiet true
        cvx_precision best
    
        variable Y(nx, nx) symmetric
        variable gammasquare

        minimize(gammasquare)

        Aa * Y + Y * Aa' + Bv * Bv' < 0
        [Y,       Y * Cc'; 
         Cc * Y,  gammasquare * eye(nz)] > 0
        Y > 0
    cvx_end
    warning('on', 'CVX:StrictInequalities');
    
    %G = Z * inv(Y);
    gval = sqrt(gammasquare);
    jval = gval / sqrt(alpha);

end
