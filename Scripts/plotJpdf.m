dat=load('plt05000/CM_temp.dat');

nVars=4;
oVar=0; % offset (counting from 0)

iT=1;
iSum=2;
iSSq=iSum+nVars;
iAvg=iSSq+nVars;
iStd=iAvg+nVars;
iN=iStd+nVars;
ip=iN+1;

grey=[0.6 0.6 0.6];

x=dat(:,iT);
a=dat(:,iAvg+oVar);
s=dat(:,iStd+oVar);

figure(1);
hold off;
fillStdDev(x,a-s,a+s,grey,0);
hold on;
plot(x,a,'k-','LineWidth',2);
xlabel('temperature (K)');
ylabel('fuel mass fraction');
