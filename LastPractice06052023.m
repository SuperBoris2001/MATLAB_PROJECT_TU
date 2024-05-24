
% строим годограф
count = 1001;
w = linspace(0,10,count);
X = zeros(1,count);
Y = zeros(1,count);

% a = [2 0 -9 0 3];
% b = [0 0 6 0 1];
% c = [2 0 1 0 1];
% d = [0 0 -6 0 5];

a = [2 0 -13 0 6.5];
b = [0 0 1.5 0 2];
c = [1 0 1.5 0 0.5];
d = [0 0 -4.5 0 8]; 
tmp1 = conv(a,b) + conv(d,c);
tmp2 = conv(a,b) - conv(d,c);
corni1 = roots(tmp1)
corni2 = roots(tmp2)
om1 = corni1(2);
om2 = corni2(2);
gam2 = -(8 - 4.5 * om2^2)/(2 + 1.5 * om2^2)
gam1 = (8 - 4.5 * om2^2)/(2 + 1.5 * om2^2)

gamma = zeros(1,count);
for i = 1:count
    X(1,i) = FindXfromW(w(i));
    Y(1,i) = FindYfromW(w(i));
    gamma(1,i) = sqrt(X(1,i)^2 + Y(1,i)^2);
end


% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);

hold on
plot(X, Y);
scatter(- gam2, - gam2);
scatter(- gam1, gam1);
plot([-gam1,gam1,gam1,-gam1,-gam1],[gam1,gam1,-gam1,-gam1,gam1])
plot([-gam2,gam2,gam2,-gam2,-gam2],[gam2,gam2,-gam2,-gam2,gam2])
%plot(x,y)
xlabel('X');
ylabel('Y');
grid()

% функции поиска x и y по w





function res = FindXfromW(w)
    res = (6.5 + 2*w^4 - 13*w^2)/(0.5 + 1.5 * w^2 + w^4);
end

function res = FindYfromW(w)
    res = (8 - 4.5 * w^2)/(2 + 1.5 * w^2);
end
