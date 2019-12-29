function f = bessel_i(n)
    ans = 1;
    n = n/2;
    
    for k=1:100
        ans = ans + ((n^k)/factorial(k))^2;
    end
    f = ans;
end