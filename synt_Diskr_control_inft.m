function synt_Diskr_control_inft
%     A = [0, 1;
%          0, 0];
%     betta = 0.1;
%     A = [0, 1;
%         -1, -2*betta];
%     Bu = [0;
%           1];
%     Bv = [0;
%           1];
    k = 1;
    m = 1;
    nu = 0.1;
    A = [0 1 0 0;
     -2*k/m -nu/m k/m 0;
     0 0 0 1;
     k/m nu/m -k/m -nu/m];
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

    T = 50;
    h = 0.05;
    tticks = (0:h:T-h); 
       
    %G = [-0.577350,-1.074570];
    %G = [-4.336729e-01,-7.525470e-01];
    %G = [0,4.442609e-01,0,-8.212825e-01];
    %2.178149e-01,3.910732e-01,-3.845442e-01,-7.826596e-01
    G = [2.178149e-01,3.910732e-01,-3.845442e-01,-7.826596e-01];
    Ac = A + Bu * G;
    count = length(tticks)     
    [~, xticks] = ode45(@(t, x)(Ac * x + Bv * exp(-0.1 * t) * sin(0.5 * t / sqrt(pi))), tticks, [0.0, 0.0, 0.0, 0.0]);
    uticks = zeros(count, 1);
    filename = 'genH2mat2_x_u_t_2_1_fourth.dat';
        fid = fopen(filename, 'wt');
        fprintf(fid, 'x,u,t\n');  % header
        fclose(fid);
    for k = 1 : count
        uticks(k) = G * xticks(k, :)';
        dlmwrite(filename, [xticks(k, 1), uticks(k), tticks(k)], ...
            'delimiter', ',', 'precision', '%.6e', '-append');
    end  
    subplot(2, 1, 1)
        plot(tticks, xticks(:, 1), 'b')
    xlabel('t')
    ylabel('x_1')
    grid on
    %title({'Непрерывная система с управлением в виде линейной обратной связи по состоянию.'})
    subplot(2, 1, 2)
        plot(tticks, uticks, 'b')
        xlabel('t')
        ylabel('u')
end