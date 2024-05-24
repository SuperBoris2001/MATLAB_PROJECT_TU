A = [0, 1, 0; 0, 1,-1; 0, 0, 0]; % решение 2 части 1 задачи из домашней работы 11.03.2023 - 19.03.2023 для другого функционала J_u  
B = [0; 0; 1];
ksi_0 = [1;1;1];
A_c_s = [];
X_s = [];
J_phi_s = [];
J_u_s = [];
i = linspace(0.01, 10, 1000);
for R = i
    Q = [1, 0 , 0; 0, 1, 0; 0, 0, R];
    X = care(A, B, Q, R);
    A_c = (A - B*1/R*B'*X);
    J_phi = ksi_0'*lyap(A_c', Q)*ksi_0;
    Q_u = X*B*R^(-2)*B'*X;
    J_u = ksi_0'*lyap(A_c', Q_u)*ksi_0;
    %A_c_s = [A_c_s, A_c];
    X_s = [X_s, X];
    J_phi_s = [J_phi_s, J_phi];
    J_u_s = [J_u_s, J_u];
end
subplot(2,2,1);
  plot(J_phi_s,J_u_s);
  grid on;
  xlabel('J_{\phi}');
  ylabel('J_u')
subplot(2,2,2)
  plot(i, J_phi_s)
  grid on
  xlabel('\rho');
  ylabel('J_{\phi}')
subplot(2,2,3)  
  plot(i, J_u_s)
  grid on;
  xlabel('\rho');
  ylabel('J_u')
