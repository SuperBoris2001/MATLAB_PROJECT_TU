alpha_min = 0;

alpha_max = 1;

alpha_count = 101;
alpha_array = linspace(alpha_min, alpha_max, alpha_count);

flag = size(alpha_array);
length_alpha_array = flag(1,2);
theta_array = zeros(2,length_alpha_array); % список theta
for i = 1:size(theta_array,2)
    theta_array(1,i) = i;
    theta_array(2,i) = i;
end

load('pqfile.mat', 'p');
load('pqfile.mat', 'q');
%save('pqfile.mat','p','q')
ax1 = subplot(1, 1, 1);
plot(theta_array(1,1:end), theta_array(2,1:end));
xlabel(ax1, '\theta_1');
ylabel(ax1, '\theta_2');