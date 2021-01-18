function [lambda,var,cum_var,err,rv,f]=re_anal(a);
%[lambda,var,cum_var,err,rv,f]=re_anal(a)
%re_anal 计算nalinowski 的re函数，rev函数及其其他统计量用于决定矩阵A中的主成份数
%lambda：a'*a的特征矢量
%var：由a'*a中特征值你和的方差
%cum_var由a'a中特征值你和的积累方差,
%err：malinowski的RE函数
%rv：malinowski简化误差特征值REV
%f：malinowski 的f-检验
nfac=min(size(a))-1;
%给结果分配矢量
lambda=zeros(nfac,1);
var=zeros(nfac,1);
cum_var=zeros(nfac,1);
err=zeros(nfac,1);
rv=zeros(nfac,1);
f=zeros(nfac,1);
pr=zeros(nfac,1);
%计算自由度
[r,c]=size(a);
y=min(r,c);
x=max(r,c);
s=svd(a',0).^2;% 计算奇异值
Trace_of_a=sum(s);%计算总平方和
lambda=s(1:nfac);
var=100.0*lambda/Trace_of_a;%计算积累方差
ssq=0.0;
resid_ssq=Trace_of_a;
for i=1:nfac  %循环计算re，cum_var,rv
    ssq=ssq+s(i);
    resid_ssq=resid_ssq-s(i);
    err(i)=sqrt(resid_ssq/((r-i)*(c-i)));
    cum_var(i)=100.0*ssq/Trace_of_a;
end;
nf=min(r,c);
for i=1:nf  %计算rev
    rv(i)=s(i)/((r-i+1)*(c-i+1));%计算简化特征值的特征矢量
end;
%计算F
for i=1:nfac
    den=sum(rv(i+1:nf));
    f(i)=(nf-i)*(rv(i)/den);
end
plot(rv)
figure
plot(err)
