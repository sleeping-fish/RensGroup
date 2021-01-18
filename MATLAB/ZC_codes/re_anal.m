function [lambda,var,cum_var,err,rv,f]=re_anal(a);
%[lambda,var,cum_var,err,rv,f]=re_anal(a)
%re_anal ����nalinowski ��re������rev������������ͳ�������ھ�������A�е����ɷ���
%lambda��a'*a������ʸ��
%var����a'*a������ֵ��͵ķ���
%cum_var��a'a������ֵ��͵Ļ��۷���,
%err��malinowski��RE����
%rv��malinowski���������ֵREV
%f��malinowski ��f-����
nfac=min(size(a))-1;
%���������ʸ��
lambda=zeros(nfac,1);
var=zeros(nfac,1);
cum_var=zeros(nfac,1);
err=zeros(nfac,1);
rv=zeros(nfac,1);
f=zeros(nfac,1);
pr=zeros(nfac,1);
%�������ɶ�
[r,c]=size(a);
y=min(r,c);
x=max(r,c);
s=svd(a',0).^2;% ��������ֵ
Trace_of_a=sum(s);%������ƽ����
lambda=s(1:nfac);
var=100.0*lambda/Trace_of_a;%������۷���
ssq=0.0;
resid_ssq=Trace_of_a;
for i=1:nfac  %ѭ������re��cum_var,rv
    ssq=ssq+s(i);
    resid_ssq=resid_ssq-s(i);
    err(i)=sqrt(resid_ssq/((r-i)*(c-i)));
    cum_var(i)=100.0*ssq/Trace_of_a;
end;
nf=min(r,c);
for i=1:nf  %����rev
    rv(i)=s(i)/((r-i+1)*(c-i+1));%���������ֵ������ʸ��
end;
%����F
for i=1:nfac
    den=sum(rv(i+1:nf));
    f(i)=(nf-i)*(rv(i)/den);
end
plot(rv)
figure
plot(err)
