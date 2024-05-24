x = [-1 0 1 2];
y = [2.4; 2.2; 0.7; -1.7];

count_x = size(x,2);
count_y = count_x;


C = zeros(4,3);
for i = 1:count_x
    C(i,1) = (x(i))^2;
    C(i,2) = x(i);
    C(i,3) = 1;
end

% 1 случай
% ABC = inv(C'*C)*C'*ySS

%2 случай
xs = [0; 0; 0];
for i = 1:count_y
    K = inv(C'*C)*C'/(i);
    xs = xs + K * (y - C * xs);
end
xs