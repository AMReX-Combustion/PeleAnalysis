function fillStdDev(x,p,m,c,t);

r=size(x,1);
a=zeros(2*r,1);
for i=1:r;
a(i)=x(i);
a(r+i)=x(1+r-i);
end;

r=size(p,1);
b=zeros(2*r,1);
for i=1:r;
b(i)=p(i);
b(r+i)=m(1+r-i);
end;

if (t==1);
h=fill(b,a,c);
else;
h=fill(a,b,c);
end;
set(h,'LineStyle','none');
