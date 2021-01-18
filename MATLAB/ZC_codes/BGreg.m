function [c,residual]=BGreg(X,Y,window)
% 该程序用于已知的背景光谱最大程度拟合数据 用数据减去拟合的背景 获得扣除背景的结果
% XY 是列光谱 
% X 是用于拟合的组分
% Y 是待拟合的一系列光谱
% c 是拟合的参数
% residual 是 扣除背景后的结果
% window 是选择 拟合的范围 [1:300] 如果是全谱 就写 ：
X=[X,ones(length(X),1)];% 带一个常数拟合
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
% 用背景数据 扣除
% plot (Fit)
% hold on 
% plot (Y)
% plot(residual)
% hold off
end
