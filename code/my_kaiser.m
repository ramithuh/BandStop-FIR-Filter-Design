function f = my_kaiser(N,alpha)
    % N should be even
    f = zeros(1,N);
    d = bessel_i(alpha);
    
    for k=0:(N-1)/2-1
        f((N-1)/2+1+k) = bessel_i(beta_(alpha,k,N))/d;
        f((N-1)/2+1-k) = bessel_i(beta_(alpha,k,N))/d;
    end
    f=f.';
end