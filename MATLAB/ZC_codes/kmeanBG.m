load result1
load s
% 方法一 用谱图直接减去背景
% BG=repmat (mclust(1,:)',1,length(sms));
% sub=sms-BG;
% save sub.mat sub
% 方法二 用谱图减去 拟合的背景
Y=sms;
X=mclust(1,:)';%输入属于背景的类数
[c,sub]=BGreg(X,Y,1:300);
save sub.mat sub c