function [y,mn,s]=autoscal(x)
% autoscal �������ֵ���Ļ����б�׼��
%[y,mn,s]=autoscal(x)
%or
%[y]=autoscal(x)
[r,c]=size(x);
mn=mean(x);
s=std(x);
y=(x-mn(ones(1,r),:))./s(ones(1,r),:);