# @Author:Xinyu_Lu
# -*- coding = 'utf-8' -*-
# @Time:2020/9/19 15:02
# @File: PCA_WPS.py
# @Software: PyCharm

import numpy as np
import matplotlib.pyplot as plt

data_O2 = np.loadtxt(r'D:\BaseProjects\Data\氧还原时间序列\O2\processed-processed-O2-OCP.txt')
data_N2 = np.loadtxt(r'D:\BaseProjects\Data\氧还原时间序列\N2\processed-processed-N2-0V.txt')
data = np.concatenate((data_O2, data_N2), axis=1)

sval_nums = 3
U, Sigma, VT = np.linalg.svd(data)
restruct = np.diag(Sigma[0:sval_nums]) @ VT[0:sval_nums, :]
restruct = (restruct - restruct.mean()) / restruct.std()  # 归一化
O2_restruct = restruct[:, 0:700]
N2_restruct = restruct[:, 700:]

# 画图1
fig = plt.figure(figsize=(6, 6))
# plt.scatter(O2_restruct[0, :], O2_restruct[1, :],color='blue')
# plt.scatter(N2_restruct[0, :], N2_restruct[1, :],color='red')
ax = fig.add_subplot(111, projection='3d')

ax.scatter(O2_restruct[0, :], O2_restruct[1, :], O2_restruct[2, :],
           marker='x', color='blue', s=40, label='class 1')
ax.scatter(N2_restruct[0, :], N2_restruct[1, :], N2_restruct[2, :],
           marker='o', color='green', s=40, label='class 2')

ax.set_xlabel('1st Component')
ax.set_ylabel('2nd Component')
plt.title('PCA Analysis')
plt.legend(['O2', 'N2'])

# 画图2
fig2 = plt.figure(figsize=(6, 6))
xList = np.arange(1, 6)
yList = Sigma[:5] / Sigma.sum() * 100
plt.plot(xList, yList)
plt.scatter(xList, yList)
plt.title('Weight Analysis')
plt.xticks(np.arange(1, 6))
plt.xlabel('No.')
plt.ylabel('Weights')
for x, y in zip(xList, yList):
    plt.text(x, y + 0.3, '%.1f' % y + '%', ha='center', va='bottom', fontsize=10.5)

plt.show()
