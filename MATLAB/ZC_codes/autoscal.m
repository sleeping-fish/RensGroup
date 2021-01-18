function [y,mn,s]=autoscal(x)
% autoscal 将矩阵均值中心化及列标准化
%[y,mn,s]=autoscal(x)
%or
%[y]=autoscal(x)
[r,c]=size(x);
mn=mean(x);
s=std(x);
y=(x-mn(ones(1,r),:))./s(ones(1,r),:);