function example
  x = linspace(-2.0,2.0,1000);
  y = linspace(-2.0,2.0,1000);
  alpha = 0.5;
  T = 1;
  [X, Y] = meshgrid(x, y);
  count = length(x);
  for i = 1 : count
      for j = 1: count
          Z(i,j) = max(OkGH1(X(i,j),Y(i,j),T)*alpha, OkGH2(X(i,j),Y(i,j),T)*(1 - alpha));
             fprintf('%d\n',Z(i,j)) 
             fprintf('%d %d\n',i,j)

      end
  end
  contour(X,Y,Z, 'ShowText', 'on');
  grid on;
end

function y = OkGH2(t1,t2,T)
    sg = t2/2;
    if t2^2 + 4 * t1 > 0
        w = sqrt((t2^2 + 4 * t1)/4);
        y = (1/(2*w))*(w - sg)*(exp((w + sg)*T) - 1) + (1/(2*w))*(w + sg)*(exp((sg - w)*T) - 1) + (sg/w)*(exp(sg + w)*T - exp(sg - w)*T);
    end
    if t2^2 + 4 * t1 < 0
     w = sqrt(-(t2^2 + 4 * t1)/4);
     fi = acos((sg^2 - w^2)/sqrt((sg^2 + w^2)^2));
     n = floor((T*w + fi)/pi);
     res = exp(-fi * sg / w);
     res = res * (w * exp(pi * sg / w) + (w * cos(fi) - sg * sin(fi)) * exp((sg*fi)/w));
     res2 = w * exp(((pi - fi) * sg) / w) * (1 + exp((pi * sg) / w));
     res2 = res2 * (1 - exp((pi * sg) * (n - 1) / w)) / (1 - exp((pi * sg) / w));
     res = res + res2;
     res3 = exp((sg / w)*(pi * n - fi));
     r = sg*(w*T + fi - pi * n)/w;
     rr = exp(r);
     
     res4 = res3 * (rr * ((sg * sin(- pi * n + T * w + fi) - w * cos (- pi * n + T * w + fi)))+ w);
     res = res + res4;
     y = res/(sg^2 + w^2);
    end
end  

function y = OkGH1(t1,t2,T)
    sg = t2/2;
    if t2^2 + 4 * t1 > 0
        w = sqrt((t2^2 + 4 * t1)/4);
        y = (1/(2*w))*(1/(w + sg))*(exp((sg + w)*T) - 1) - (1/(2*w))*(1/(sg - w))*(exp((sg - w)*T) - 1);
    end
    if t2^2 + 4 * t1 < 0
     w = sqrt(-(t2^2 + 4 * t1)/4);
     n = floor(T*w/pi); 
     %fprintf('%d %d \n',t1,t2);
     res = (1 - exp(pi * sg * n / w))*(1 + exp(pi * sg / w)) / (1 - exp(pi * sg / w));
     res = res + exp(pi * sg * n / w)*(w - exp(sg*(T - pi * n / w)) * (-1)^n * (w*cos(T*w)-sg*sin(T*w)))/w;
     res = res/(sg^2 + w^2);
     y = res;
    end
end  

% function y = GH(t1,t2)
%   T = 1;
%   A = [0 1; 
%       t1 t2];
%   B = [0; 1];
%   C = [t1 t2];
%   y = integral(@(tau)(abs( C * expm(A * tau) * B)), 0, T,'ArrayValued',true);
% end




 


