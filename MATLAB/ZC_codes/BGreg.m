function [c,residual]=BGreg(X,Y,window)
% �ó���������֪�ı����������̶�������� �����ݼ�ȥ��ϵı��� ��ÿ۳������Ľ��
% XY ���й��� 
% X ��������ϵ����
% Y �Ǵ���ϵ�һϵ�й���
% c ����ϵĲ���
% residual �� �۳�������Ľ��
% window ��ѡ�� ��ϵķ�Χ [1:300] �����ȫ�� ��д ��
X=[X,ones(length(X),1)];% ��һ���������
[ry,cy]=size(Y);
[rx,cx]=size(X);
residual=zeros(ry,cy);
c=zeros(cx,cy);
x=X(window,:);
y=Y(window,:);
for ii=1:cy
[c(:,ii)]=lsqnonneg(x,y(:,ii));
Fit=X*c(:,ii);
residual(:,ii)=Y(:,ii)-Fit;
% �ñ������� �۳�
% plot (Fit)
% hold on 
% plot (Y)
% plot(residual)
% hold off
end
