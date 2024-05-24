A = [0, 1; 0, 0];
Bu = [0; 1];
Bv = [0; 1];
C1 = [1, 0];
alphasquare = 0.5;
THETA = [1,1];
x = [1,1];
I = eye(2);
h = 0.01;
cvx_begin
 variable Y(2,2, 100) symmetric;
 variable gammasquare;
 minimize(gammasquare);
 subject to
 [Y(:,:,1) , (C1 * Y(:,:,1))'; C1 * Y(:,:,1), alphasquare*gammasquare] >= 0;
 [Y(:,:,1), (THETA * Y(:,:,1))'; THETA*Y(:,:,1), (1 - alphasquare)^2 * gammasquare] >= 0;
 for i = 2:100
 [(Y(:,:,i)-Y(:,:,i-1))/h + (A + Bu*THETA)*Y(:,:,i-1) + Y(:,:,i-1)*(A + Bu*THETA)', Bv; Bv', - 1] <= 0;
 [Y(:,:,i) , (C1 * Y(:,:,i))'; C1 * Y(:,:,i), alphasquare*gammasquare] >= 0;
 [Y(:,:,i), (THETA * Y(:,:,i))'; THETA*Y(:,:,i), (1 - alphasquare)^2 * gammasquare] >= 0;
 end
cvx_end

%load('levitation_system_01.mat',"tticks","x1ticks","x2ticks","x3ticks","uticks");
% load(levitation_system_01,"-mat");