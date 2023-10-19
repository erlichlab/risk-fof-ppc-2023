function o = mseutility(b, x, y)

 u = @(b,x)(b(2).*x.^b(1));
 b0 = b;

 if b(1) < 0
     b(1) =0;
 end
 if b(2) < 0
     b(2) =0;
 end

 ey = u(b,x);

 o = sum((y - ey).^2) + sum(abs(b0)); % L1 penalty on parameters

