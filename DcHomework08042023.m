Source = load('dc_motor.mat');
uticks = Source.uticks;
tticks = Source.tticks;
xticks = Source.xticks;
flag = length(tticks);
LIST_C = cell(1,flag);
LIST_T = cell(1,flag);
LIST_Y = cell(1,flag);
LIST_X = cell(1,flag);
h = 0.1;
SUM_C = zeros(5,5);
LIST_X{1} = [7.5;7.5;7.5;7.5;7.5];

for i = 2:length(tticks)
    LIST_C{i-1} = [xticks(i,1)-xticks(i-1,1), - xticks(i,2),0 ,0 ,0 ;
        0 ,0 ,xticks(i,2)-xticks(i-1,2), xticks(i-1,2)*h ,xticks(i-1,1)*h]; 
    SUM_C = SUM_C + LIST_C{i-1}'*LIST_C{i-1};
    LIST_T{i-1} = 2*LIST_C{i-1}'*LIST_C{i-1} + SUM_C; 
    LIST_Y{i-1} = [0;uticks(i-1)*h];
    LIST_X{i} = inv(LIST_T{i-1});
    LIST_X{i} = LIST_X{i} * LIST_C{i-1}'*(LIST_Y{i-1} - LIST_C{i-1}*LIST_X{i-1});
    LIST_X{i} = LIST_X{i-1} - LIST_X{i};
end


ax1 = subplot(4,2,1);
plot(ax1,tticks, xticks);
xlabel(ax1,'t');
ylabel(ax1,'x');


ax2 = subplot(4,2,2);
plot(ax2,tticks, uticks);
xlabel(ax2,'t');
ylabel(ax2,'u');


ax3 = subplot(4,2,3);
J = zeros(1,flag);
for i = 1:flag
    J(i) = LIST_X{i}(1); 
end
plot(ax3,tticks, J);
xlabel(ax3,'t');
ylabel(ax3,'J');


ax4 = subplot(4,2,4);
Kt = zeros(1,flag);
for i = 1:flag
    Kt(i) = LIST_X{i}(2); 
end
plot(ax4,tticks, Kt);
xlabel(ax4,'t');
ylabel(ax4,'Kt');

ax5 = subplot(4,2,5);
L = zeros(1,flag);
for i = 1:flag
    L(i) = LIST_X{i}(3); 
end
plot(ax5,tticks, L);
xlabel(ax5,'t');
ylabel(ax5,'L');

ax6 = subplot(4,2,6);
R = zeros(1,flag);
for i = 1:flag
    R(i) = LIST_X{i}(4); 
end
plot(ax6,tticks, R);
xlabel(ax6,'t');
ylabel(ax6,'R');

ax7 = subplot(4,2,7);
Km = zeros(1,flag);
for i = 1:flag
    Km(i) = LIST_X{i}(5); 
end
plot(ax7,tticks, Km);
xlabel(ax7,'t');
ylabel(ax7,'Km');


