function genH4_inf_res

%     задаем матрицы системы
%     A = [0, 1;
%          0, 0];
    
%     betta = 0.1;
%     A = [0, 1;
%     -1, -2*betta];
%     Bu = [0;
%           1];
%     Bv = [0;
%           1];
%     C1 = [1, 0];
%     D1 = 0;
%     C2 = [0, 0];
%     D2 = 1;  

    k = 1;
    m = 1;
    nu = 0.1;
%     A = [0 1 0 0;
%      -2*k/m -nu/m k/m 0;
%      0 0 0 1;
%      k/m nu/m -k/m -nu/m];    
    A = [0 1 0 0;
     -2*k/m 0 k/m 0;
     0 0 0 1;
     k/m 0 -k/m 0];
    Bu = [0;
          0;
          0;
          1];
    Bv = [0;
          1;
          0;
          0];
    C1 = [-1, 0, 1, 0];
    D1 = 0;
    C2 = [0, 0, 0, 0];
    D2 = 1;
    
    lamlist = 0.02 : 0.02 : 0.98;
    count = length(lamlist);
    vals_list = zeros(count, 7);
    k = 1;
    
    filename = 'genH2mat4inf_res_3.dat';
        fid = fopen(filename, 'wt');
        fprintf(fid, 'alpha_l,J1,J2,G1,G2\n');  % header
        fclose(fid);
    
    for lam = lamlist
        %
        [~, G, gval] = genH2(A, Bu, Bv, C1, D1, C2, D2, lam);
        % формируем матрицы замкнутой системы
        Ac = A + Bu * G; 
        vals_list(k, 1) = lam;
        vals_list(k, 2) = gen_H2_norm(Ac, Bv, C1 + D1 * G);
        vals_list(k, 3) = gen_H2_norm(Ac, Bv, C2 + D2 * G);
        vals_list(k, 4) = G(1,1);
        vals_list(k, 5) = G(1,2);
        vals_list(k, 6) = G(1,3);
        vals_list(k, 7) = G(1,4);
        fprintf('--------------------------------------------------------------------------------------------------------\n');
        fprintf('alpha is %4.6f | J1 is %4.6f |J2 is %4.6f|G1 is %4.6f |G2 is %4.6f |G3 is %4.6f |G4 is %4.6f',vals_list(k,1),vals_list(k,2),vals_list(k,3),vals_list(k,4),vals_list(k,5),vals_list(k,6),vals_list(k,7));
        fprintf('\n');
        % ;val_list(k,4);val_list(k,5)   
        dlmwrite(filename, [vals_list(k,1);vals_list(k,2);vals_list(k,3);vals_list(k,4);vals_list(k,5);vals_list(k,6);vals_list(k,7)]', ...
            'delimiter', ',', 'precision', '%.6e', '-append');
       
        
        k = k + 1;
    end
    
    plot(vals_list(:, 2), vals_list(:, 3))
    
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

function gval = gen_H2_norm(Ac, Bv, Cc)
   
    % определяем размер
    nx = size(Ac, 1);
    nz = size(Cc, 1);
       
    warning('off', 'CVX:StrictInequalities');
    cvx_begin sdp
        cvx_solver sdpt3
        cvx_precision best
        cvx_quiet true
        
        variable Y(nx, nx) symmetric
        variable gammasquare

        minimize(gammasquare)

        Ac * Y + Y * Ac' + Bv * Bv' < 0
        [Y,       Y * Cc'; 
         Cc * Y,  gammasquare * eye(nz)] > 0
        Y > 0
    cvx_end
    warning('on', 'CVX:StrictInequalities');
    
    gval = sqrt(gammasquare);

end
