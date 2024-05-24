function y = OkGH(t1,t2)
   sg = t2/2;
   w = sqrt(-(t2^2 + 4 * t1^2)/4);
   fi = acos((sg^2 - w^2)/sqrt(sg^2 + w^2));
   n = floor(T/w);
   res = exp^(-fi * sg / w);
   res = res * (w * exp ^ (pi * sg / w) + (w * cos(fi) - sg * sin(fi)) * exp ((1 + sqrt(5))/(2 * w)));
   res2 = w * exp ^ ((pi * sg) / w - fi) * (1 + exp ^ ((pi * sg) / w));
   res2 = res2 * (1 - exp ^ ((pi * sg) * (n - 1) / w)) / (1 - exp ^ ((pi * sg) / w));
   res = res + res2;
   res3 = exp ^ ((sg / w)*(pi * n - fi));
   res3 = res3 * (exp^(sg * (-2 * pi * n + 2 * T * w + 1 + sqrt(5))/(2 * w)) * (sg * sin(- pi * n + T * w + fi) - w * cos (- pi * n + T * w + fi)) + w);
   res = res + res3;
   y = res;
end  
    
 