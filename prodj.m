function prodj()
 x0 = [2.5,3.08];
 [x, f] = fminsearch(@mixOkGH,x0);
 fprintf('%d %d \n',x,f);
end


function y = mixOkGH(t)
   T = 5;
 
   alpha = 0.15;
   y = max(OkGH1(t, T) * alpha, OkGH2(t, T) * (1 - alpha));
   
end

function y = OkGH1(t,T)%C_1 = [1, 0]
   t1 = t(1);
   t2 = t(2);
   sg = t2/2;
   if t2^2 + 4 * t1 > 0
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
   if t2^2 + 4 * t1 > 0
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

