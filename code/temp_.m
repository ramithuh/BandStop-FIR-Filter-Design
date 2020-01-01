
%%Function to Calculate the kaiser window
function f = my_kaiser(N,alpha)
    % N should be odd
    f = zeros(1,N);
    d = bessel_i(alpha);
    
    for k=0:(N-1)/2-1
        f((N-1)/2+1+k) = bessel_i(beta_(alpha,k,N))/d;
        f((N-1)/2+1-k) = bessel_i(beta_(alpha,k,N))/d;
    end
    f=f.';
end

%%Function to calculate the modified Bessel function of the first kind
function f = bessel_i(n)
    ans = 1;
    n = n/2;
    
    for k=1:100
        ans = ans + ((n^k)/factorial(k))^2;
    end
    f = ans;
end

%%Function to calculate beta value for a given alpha
function f = beta_(alpha,n,N)
    f = alpha*(1-((2*n)/(N-1))^2)^0.5;
end