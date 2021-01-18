load s
s=sms;
nk=7;%设定聚类的类数
[r,c]=size(s);
s=s';
opts=statset('Display','final');
[idx,ctrs]=kmeans(s,nk,'Replicates',5,'Options',opts);%可以用idx数值重构图像
for i=1:nk
clust=s(idx==i,:);% 找到所有的属于i类的谱图
mclust(i,:)=mean(clust);%获得每一类的平均谱图，一行为一张谱图。
end
save('result1.mat','mclust','idx','ctrs')
H=400;%水平pixel
V=101;%竖直pixel
pidx=reshape(idx(:),H,V);
pidx=pidx';
figure
pcolor (pidx)
shading flat % 无轮廓线 无差值
axis equal tight % xy 单位长度一致
axis ij % y轴反向
figure
plot(mclust')
    