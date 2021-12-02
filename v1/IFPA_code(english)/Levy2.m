function L = Levy2(d,iter,ger)
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;%Standard normal distribution of u
    v=randn(1,d);%Standard normal distribution of v
    step=u./abs(v).^(1/beta);
    a = (1 - (iter./ger).^2).^0.5;
    L = 0.01*a*step;
end