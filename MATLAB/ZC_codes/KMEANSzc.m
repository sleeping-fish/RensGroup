load s
s=sms;
nk=7;%�趨���������
[r,c]=size(s);
s=s';
opts=statset('Display','final');
[idx,ctrs]=kmeans(s,nk,'Replicates',5,'Options',opts);%������idx��ֵ�ع�ͼ��
for i=1:nk
clust=s(idx==i,:);% �ҵ����е�����i�����ͼ
mclust(i,:)=mean(clust);%���ÿһ���ƽ����ͼ��һ��Ϊһ����ͼ��
end
save('result1.mat','mclust','idx','ctrs')
H=400;%ˮƽpixel
V=101;%��ֱpixel
pidx=reshape(idx(:),H,V);
pidx=pidx';
figure
pcolor (pidx)
shading flat % �������� �޲�ֵ
axis equal tight % xy ��λ����һ��
axis ij % y�ᷴ��
figure
plot(mclust')
    