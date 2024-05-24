v = linspace(1,10,1000);
w = linspace(1,10,1000);
tticks = linspace(1,10,1000);
flag = length(tticks) - 1;
xticks = linspace(1,10,1000);
yticks = linspace(1,10,1000);
A = 1;
xticks(1) = rand;
yticks(1) = rand;
for i = 2:flag
    yticks(i) = xticks(i - 1) + normrnd(0, 0.01);
    xticks(i) = xticks(i - 1) + normrnd(0, 0.01);    
end

ax1 = subplot(3,1,1);
plot(tticks, yticks);
xlabel(ax1,'t');
ylabel(ax1,'y');
subplot(3,1,2);

ax2 = subplot(3,1,2);
plot(tticks, xticks);
xlabel(ax2,'t');
ylabel(ax2,'x');

ax3 = subplot(3,1,3);
plot(xticks, yticks);
xlabel(ax3,'x');
ylabel(ax3,'y');