load result1
load s
% ����һ ����ͼֱ�Ӽ�ȥ����
% BG=repmat (mclust(1,:)',1,length(sms));
% sub=sms-BG;
% save sub.mat sub
% ������ ����ͼ��ȥ ��ϵı���
Y=sms;
X=mclust(1,:)';%�������ڱ���������
[c,sub]=BGreg(X,Y,1:300);
save sub.mat sub c