function da1_c=diff_center(a1_c)
N=length(a1_c);
h=1/(N-1);

a1_ce=[2*a1_c(1)-a1_c(2);a1_c;2*a1_c(end)-a1_c(end-1)];
da1_c=(a1_ce(3:end)-a1_ce(1:end-2))/(2*h);
