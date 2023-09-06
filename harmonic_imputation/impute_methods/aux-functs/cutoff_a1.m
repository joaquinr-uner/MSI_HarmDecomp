function y=cutoff_a1(x,a)
y=(-cos((x-a)*pi/a)+1)/2;
b=(sign(x)+1)/2;
y=b.*y+(1-b);
c=(sign(a-x)+1)/2;
y=y.*c;
