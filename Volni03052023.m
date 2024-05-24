syms w h t s
f = (exp(-(w - h)^2/2 - (w + h)^2/8 - 1i * w * t + 1i * h * s));
res = int(int(f,h,-inf, inf),w, -inf, inf) 
