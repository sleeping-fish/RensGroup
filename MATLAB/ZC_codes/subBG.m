clc
clear
load DataWN.txt
load BG.txt
[r,c]=size(DataWN);
bg=BG(:,2);
bg=repmat(bg,1,c-1);
WN=DataWN(:,1);
sbg=DataWN(:,2:end)-bg;
sDataWN=[WN,sbg];%�ڵ�һ�м���ǿ����Ϣ
dlmwrite ('sDataWN.txt', sDataWN,'\t')