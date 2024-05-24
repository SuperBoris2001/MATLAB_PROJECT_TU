function synt_Diskr_control
    Source = load('genH2mat2_theta_second_test.dat','-ascii');
    %alpha,J1,J2,theta_1,theta_2
    theta_1 = Source(:,4);
    theta_2 = Source(:,5);
    A = [0, 1;
         0, 0];
    % betta = 0.1;
    % A = [0, 1;
    %     -1, -2*betta];
    Bu = [0;
          1];
    Bv = [0;
          1];
    old_tticks = Source(:,6);
    h = old_tticks(2) - old_tticks(1);
    old_count = length(old_tticks);
    xticks = [];
    tticks = [];
    uticks = [];
    for k = 1:old_count
        G = [theta_1(k),theta_2(k)];
        low_lim = old_tticks(k);
        high_lim = old_tticks(k) + h;
        tiks = (low_lim:0.001:high_lim);
        Ac = A + Bu * [theta_1(k),theta_2(k)];
        if k == 1
           [tk1, xk1] = ode45(@(t, x)(Ac * x + Bv * exp(-0.1 * t) * sin(0.5 * t / sqrt(pi))),tiks, [0, 0]); 
        else
           [tk1, xk1] = ode45(@(t, x)(Ac * x + Bv * exp(-0.1 * t) * sin(0.5 * t / sqrt(pi))),tiks, xticks(end,:));
        end
        xticks = [xticks;xk1];
        tticks = [tticks;tk1];
        res = xk1*G';
        uticks = [uticks;res];
    end
    count = length(tticks);
    filename = 'genH2mat2_x_u_t_1_second.dat';
        fid = fopen(filename, 'wt');
        fprintf(fid, 'x,u,t\n');  % header
        fclose(fid);
    dlmwrite(filename, [xticks(:,1) uticks tticks], ...
         'delimiter', ',', 'precision', '%.6e', '-append');
    subplot(2, 1, 1)
    plot(tticks, xticks(:,1), 'b')
    xlabel('t')
    ylabel('x_1')
    grid on
    subplot(2, 1, 2)
    plot(tticks, uticks, 'b')
    xlabel('t')
    ylabel('u') 
end

