function example2
    % задаем кол-во точек в сетке
    count = 5000;
    x = linspace(-5.0, 5.0, count + 1);
    y = linspace(-5.0, 5.0, count + 1);
   
    % задаем параметр свертки
    alpha = 0.5;
    T = 10;
    [X, Y] = meshgrid(x, y);
    Z = zeros(count + 1);
    for i = 1 : count + 1
        for j = 1 : count + 1
            Z(i,j) = max(OkGH1([X(i, j),Y(i, j)], T) * alpha, ...
                         OkGH2([X(i, j), Y(i, j)], T) * (1 - alpha));
            %fprintf('%d \n',Z(i,j));       
        end
    end
    pnts_lst = zeros(2, count + 1);
    for k = 1 : count + 1
        pnts_lst(2, k) = y(k);
        pnts_lst(1, k) = -0.25 * y(k)^2;
    end
        
    contour(X, Y, Z, [0.1, 0.2, 0.3, 0.4, 1.5], 'ShowText', 'on');
    hold on;
    plot(pnts_lst(1, :), pnts_lst(2, :), 'k','LineWidth', 2.0);
    grid on;

end
function y = OkGH1(t,T)%C_1 = [1, 0]
   t1 = t(1);
   t2 = t(2);
   sg = t2/2;
   if t2^2 + 4 * t1 >= 0
       w = sqrt((t2^2 + 4 * t1)/4);
       y = abs((1/(2*w))*(1/(w + sg))*(exp((sg + w)*T) - 1) - (1/(2*w))*(1/(sg - w))*(exp((sg - w)*T) - 1));
       
   elseif t2^2 + 4 * t1 < 0
       w = sqrt(-(t2^2 + 4 * t1)/4);
       n = floor(T*w/pi);
       
       res = (1 - exp(pi * sg * n / w))*(1 + exp(pi * sg / w)) / (1 - exp(pi * sg / w));
       res = res + exp(pi * sg * n / w)*(w - exp(sg*(T - pi * n / w)) * (-1)^n * (w*cos(T*w)-sg*sin(T*w)))/w;
       res = res/(sg^2 + w^2);
       
       y = res; 
   end
end

function y = OkGH2(t,T)%C_2 - [t1, t2]
   t1 = t(1);
   t2 = t(2);
   sg = t2/2;
   if t2^2 + 4 * t1 >= 0
       w = sqrt((t2^2 + 4 * t1)/4);
       y = abs((1/(2*w))*(w + sg)*(exp((w + sg)*T) - 1) + (1/(2*w))*(w - sg)*(exp((sg - w)*T) - 1));
       if y < 0
        fprintf('%d\n',y); 
       end
   elseif t2^2 + 4 * t1 < 0
       w = sqrt(-(t2^2 + 4 * t1)/4);
       fi = acos((sg^2 - w^2)/(sg^2 + w^2));
       n = floor((T*w + fi)/pi);
       res1 = exp((pi-fi) * sg / w);
       res1 = res1 * (1 + (1 + exp(pi*sg/w))*(1 - exp(pi*sg*(n-1)/w))/(1 - exp(pi*sg/w)));
       res2 = exp(sg*(pi-fi)/w)/w;
       res2 = res2 * (exp(sg*(w*T + fi - pi*n)/w)*(sg*sin(w*T + fi - pi*n) - w* cos(w*T + fi - pi*n)) + w);
       res3 = (w*cos(fi) - sg*sin(fi));
       res = res1 + res2 + res3;       
       
       y = res;
   end
end  
% function y = OkGH2(t1,t2,T)
%    sg = t2/2;
%    if t2^2 + 4 * t1 >= 0
%        w = sqrt((t2^2 + 4 * t1)/4);
%        y = (1/(2*w))*(w - sg)*(exp((w + sg)*T) - 1) + (1/(2*w))*(w + sg)*(exp((sg - w)*T) - 1) + (sg/w)*(exp(sg + w)*T - exp(sg - w)*T);
%    else
%        w = sqrt(-(t2^2 + 4 * t1)/4);
%        fi = acos((sg^2 - w^2)/sqrt((sg^2 + w^2)^2));
%        n = floor((T*w + fi)/pi);
%        res = exp(-fi * sg / w);
%        res = res * (w * exp(pi * sg / w) + (w * cos(fi) - sg * sin(fi)) * exp((sg*fi)/w));
%        res2 = w * exp(((pi - fi) * sg) / w) * (1 + exp((pi * sg) / w));
%        res2 = res2 * (1 - exp((pi * sg) * (n - 1) / w)) / (1 - exp((pi * sg) / w));
%        res = res + res2;
%        res3 = exp((sg / w)*(pi * n - fi));
%        r = sg*(w*T + fi - pi * n)/w;
%        rr = exp(r);
%    
%        res4 = res3 * (rr * ((sg * sin(- pi * n + T * w + fi) - w * cos (- pi * n + T * w + fi)))+ w);
%        res = res + res4;
%        y = res/(sg^2 + w^2);
%    end
% end  
% 
% function y = OkGH1(t1,t2,T)
%    sg = t2/2;
%    if t2^2 + 4 * t1 >= 0
%        w = sqrt((t2^2 + 4 * t1)/4);
%        y = (1/(2*w))*(1/(w + sg))*(exp((sg + w)*T) - 1) - (1/(2*w))*(1/(sg - w))*(exp((sg - w)*T) - 1);
%    else
%        w = sqrt(-(t2^2 + 4 * t1)/4);
%        n = floor(T*w/pi);
%        %fprintf('%d %d \n',t1,t2);
%        res = (1 - exp(pi * sg * n / w))*(1 + exp(pi * sg / w)) / (1 - exp(pi * sg / w));
%        res = res + exp(pi * sg * n / w)*(w - exp(sg*(T - pi * n / w)) * (-1)^n * (w*cos(T*w)-sg*sin(T*w)))/w;
%        res = res/(sg^2 + w^2);
%        y = res;
%    end
% end

