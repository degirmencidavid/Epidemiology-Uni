function y=SEIRDfunction(t,x) %¼ÆËãº¯Êý
    N = 10000; %×ÜÈËÊý
    y = [ -5.*x(1).*x(3)./N...%S
        5.*x(1).*x(3)./N-(1/10).*x(2)...%E
        (1/10).*x(2)-(1/3).*x(3)...%I
        0.99*(1/3).*x(3)...%R
        (0.01)*(1/3).*x(3)]';%D
    %¹¹ÔìÓÃÓÚ½âÎ¢·Ö·½³ÌµÄÁÐÏòÁ¿
end