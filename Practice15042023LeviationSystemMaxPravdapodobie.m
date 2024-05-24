Source = load('levitation_system_01');
tticks = Source.tticks;
uticks = Source.uticks;
xticks = Source.xticks;
a = 7.5;
A = [0, 1, 0; 1, 0, 1; 0, -1, -a];
flag = length(tticks);

Sigma0 = xticks(1,1:3);

LIST_Sigma_minus = cell(1,flag);
LIST_Sigma_plus = cell(1,flag);
LIST_X_minus = cell(1,flag);
LIST_X_plus = cell(1, flag);

    
%
% B = [0; 0; 1];
% C = [0, 0, 1];

ax1 = subplot(4,2,1);
plot(ax1,tticks, xticks);
xlabel(ax1,'t');
ylabel(ax1,'x');


ax2 = subplot(4,2,2);
plot(ax2,tticks, uticks);
xlabel(ax2,'t');
ylabel(ax2,'u');
