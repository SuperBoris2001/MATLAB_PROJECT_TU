A = [0 0 1 0;
    0 0 0 1;
    74.49823574469526 -38.256590310670866 -56.04604080240516 86.34286705234221;
    -114.76977093201259 101.34592586957349 86.34286705234221 -228.73177490708963];
B = [0;
    0;
    112.0919;
    -172.6856];

% g = 9.81;
% m1 = 0.127;
% m2 = m1;
% I1 = 1.2 * 10^(-3);
% I2 = I1;
% l1 = 0.178;
% l2 = l1;
% L1 = 0.356;
% L2 = L1;
% K1 = 1.4 * 10^(-3);
% K2 = K1;
% Kf = 8.48 * 10^(-3);
% Ks = 0.028;
% A11 = I1 + m1 * l1^2 + m2 * L1^2;
% A12 = m2 * l2 * L1;
% A22 = I2 + m2 * l2^2;
% A = [0 0 1 0; 
%     0 0 0 1;
%     (A22 * (m1 * l1 + m2 * L1) * g)/(A11 * A22 - A12^2) -(A12 * m2 * g * l2)/(A11 * A22 - A12^2) -(A22 * (Kf * Ks + K1))/(A11 * A22 - A12^2) (A12 * K2)/(A11 * A22 - A12^2);
%     -(A12 * (m1 * l1 + m2 * L1) * g)/(A11 * A22 - A12^2) (1 + A12^2/(A11 * A22 - A12^2))*(m2 * g * l2)/(A22) (A12 * (Kf * Ks + K1))/(A11 * A22 - A12^2) -(1 + A12^2/(A11 * A22 - A12^2))* K2 / A22];
% 
% B = [0;
%     0;
%     (A22 * Kf)/(A11 * A22 - A12^2);
%     -(A12 * Kf)/(A11 * A22 - A12^2)];







C = [1 1 0 0];

% -265.06819578, -21.41386813, -1+2i, -1-2i
% -18.85486976, -4.5500792, 4.46864423, 12.11231171
K = place(A', C', [-265.06819578, -21.41386813, -1,-2]);
L = K'

H = -place(A, B, [-265.06819578, -21.41386813, -1, -2])
%РЈРїСЂР°РІР»РµРЅРёРµ РїРѕ СЃРѕСЃС‚РѕСЏРЅРёСЋ РЅР°Р±Р»СЋРґР°С‚РµР»СЏ
G = [A, B*H; L*C, A + B * H - L * C]

%РќР°Р±Р»СЋРґР°С‚РµР»СЊ
%G = [A, zeros(4); L*C, A - L * C];
%[ttiks, xticks] = ode45(@(t,x)(G*x), (0: 0.0001: 10), [0.01, 0.01, 0, 0, 0, 0, 0, 0]);
[ttiks, xticks] = ode45(@(t,x)(G*x), (0: 0.0001: 10), [0.6, 0.1, 0, 0, 0, 0, 0, 0]);

% линейная система
subplot(6, 2, 1)
plot(ttiks,xticks(:, 1),'b', ttiks,xticks(:, 5),'r')
xlabel('t')
ylabel('\theta_1')
subplot(6, 2, 2)
plot(ttiks,xticks(:, 2),'b', ttiks,xticks(:, 6),'r')
xlabel('t')
ylabel('\theta_2')
subplot(6, 2, 3)
plot(ttiks, xticks(:, 3),'b', ttiks,xticks(:, 7),'r')
xlabel('t')
ylabel({'$\dot{\theta}_1$'},'Interpreter','latex')
subplot(6, 2, 4)
plot(ttiks, xticks(:, 4),'b', ttiks,xticks(:, 8),'r')
xlabel('t')
ylabel({'$\dot{\theta}_2$'},'Interpreter','latex')


count = length(ttiks);
zticks = zeros(count, 1);
uticks = zeros(count, 1);
for k = 1 : count
zticks(k) =C * (xticks(k,1:4) - xticks(k,5:8))';
uticks(k) = H * (xticks(k,5:8))';
end

% РќРµРІСЏР·РєР° (РґР»СЏ РЅР°Р±Р»СЋРґР°С‚РµР»СЏ)
subplot(6,2,5)
plot(ttiks, zticks)
xlabel('t')
ylabel({'$\varepsilon$'},'Interpreter','latex')

% РЈРїСЂР°РІР»РµРЅРёРµ РїРѕ СЃРѕСЃС‚РѕСЏРЅРёСЋ РЅР°Р±Р»СЋРґР°С‚РµР»СЏ
subplot(6,2,6)
plot(ttiks, uticks)
xlabel('t')
ylabel('u')

%РџСЂРѕРІРµСЂРєР° СЃРѕР±СЃС‚РІРµРЅРЅС‹С… С‡РёСЃРµР»
eig(A-L*C)
eig(A+B*H)



% нелинейная система

% синтезируем ассимптотический наблюдатель для нелинейной системы
% [ttiks, x_n] = ode45(@(t,X)to_ode_1(t, X, A, B, H, L, C, Kf, Ks, K1, K2, A11, A12, A22, g, l1, l2, L1, L2, m1, m2), (0: 0.1: 10), [0.6, 0.1, 0, 0, 0, 0, 0, 0]);
% 
% subplot(6,2,6)
% plot(ttiks, x_n(:,1),  'b', ttiks, x_n(:, 5), 'r')
% xlabel('t')
% ylabel({'$\theta_1$'},'Interpreter','latex')
% 
% subplot(6,2,7)
% plot(ttiks, x_n(:,2),  'b', ttiks, x_n(:, 6), 'r')
% xlabel('t')
% ylabel({'$\theta_2$'},'Interpreter','latex')
% 
% subplot(6,2,8)
% plot(ttiks, x_n(:,3),  'b', ttiks, x_n(:, 7), 'r')
% xlabel('t')
% ylabel({'$\dot{\theta_1}$'},'Interpreter','latex')
% 
% subplot(6,2,9)
% plot(ttiks, x_n(:,4),  'b', ttiks, x_n(:, 8), 'r')
% xlabel('t')
% ylabel({'$\dot{\theta_2}$'},'Interpreter','latex')
% 
% function dxdt = to_ode_1(t, X, A, B, H, L, C, Kf, Ks, K1, K2, A11, A12, A22, g, l1, l2, L1, L2, m1, m2)
%     ksi_1 = X(5);
%     ksi_2 = X(6);
%     ksi_3 = X(7);
%     ksi_4 = X(8);
%     ksi = [ksi_1; ksi_2; ksi_3; ksi_4];
%     dksidt = A*ksi + B*H*ksi + L*C*(X(1:4)-ksi);
%     dksi_1dt = dksidt(1);
%     dksi_2dt = dksidt(2);
%     dksi_3dt = dksidt(3);
%     dksi_4dt = dksidt(4);
%     u = H*ksi;
%     
%     dx_1dt = X(3);
%     dx_2dt = X(4);  
% %     dx_3dt = (-A12*dx_2dt^2 * sin(X(1)-X(2))+g*sin(X(1))*(m1*l1+m2*L1)-A12^2 * dx_1dt^2 * sin(2*(X(1)-X(2)))/(2*A22)-A12*cos(X(1)-X(2))*m2*g*l2*sin(X(2))/A22 + Kf*(u-Ks*dx_1dt)-K1*dx_1dt+A12*cos(X(1)-X(2))*K2*dx_2dt/A22) / (A11 - A12^2 * cos(X(1)-X(2))^2 / A22);
% %     dx_4dt = (-A12*dx_3dt*cos(X(1)-X(2)) + A12*dx_1dt^2 *sin(X(1)-X(2))+m2*g*l2*sin(X(2))-K2*dx_2dt) / A22;
%     
%     M = Kf * (u - Ks * X(3));
%     tmp2 = (A11 * A22 - A12^2 * (cos(X(1) - X(2)))^2);
%     tmp1 = (M - K1 * X(3) + A12 * K2 * X(4) * cos(X(1) - X(2)) / A22 - A12 * (X(4))^2 * sin(X(1) - X(2)) + g * sin(X(1)) * (m1 * l1 + m2 * L1) - A12 * m2 * g * l2 * sin(X(2)) * cos(X(1) - X(2))/A22 );
%     %   
%     dx_3dt = tmp1 / tmp2;
%     tmp3 = (- K2 * cos(X(1) - X(2)) * X(4) / A22 + m2 * g * l2 * sin(X(2)) * cos(X(1) - X(2)) / A22 + A12 * (X(3))^2 * sin(X(1) - X(2)) * cos(X(1) - X(2)) / A12 - A12 * X(1) * (cos(X(1) - X(2))^2 / A22));
%     tmp4 = (cos(X(1) - X(2)));
%     dx_4dt = tmp3/tmp4;
%     dxdt = [dx_1dt; dx_2dt; dx_3dt; dx_4dt; dksi_1dt; dksi_2dt; dksi_3dt; dksi_4dt];
% end

