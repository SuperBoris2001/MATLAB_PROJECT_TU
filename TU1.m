A = [-1, 1; 0, -2];
C = [1, 0];
L = [-2.8; 4.61];
A=[A, zeros(2);L*C,A-L*C];

eig(A)

[tticks, xticks] = ode45(@(t, x)(A*x), [0,30],[1, 1,0,0]);
subplot(3,1,1);
plot(tticks,xticks(:,1),'b',tticks,xticks(:,3),'r');
grid;
subplot(3,1,2);
plot(tticks,xticks(:,2),'b',tticks,xticks(:,4),'r');
grid;
count=length(tticks);
zticks=zeros(count,1);
for k=1:count
zticks(k)=C*(xticks(k,1:2)-xticks(k,3:4))';
end
subplot(3,1,3);
plot(tticks, zticks);
grid;
A = [-1, 1; 0, -2];
C = [1, 1; 0, 1];

roots = [-2 + 1i, -2 - 1i];
hi_star = poly(roots);
hi_znach = polyval(hi_star, A');

k = - hi[0,1]*inv(C)*hi_znach;
%k= place(A',C',[-0.1+1i;-0.1-1i]);
k  