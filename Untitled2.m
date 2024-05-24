a = 7.5;
function HomeProject()
 fprintf('%d %d \n',x,f);
end
function y = matrr(x,u)
  y = [x(2), 0.5*((1 + x(3))^2/(1 - x(1))^2 - 1),-(1 + x(3))/(1 - x(1))*x(2) - a*(1-x(1))*x(3) + (1 - x(1))*u];
end
