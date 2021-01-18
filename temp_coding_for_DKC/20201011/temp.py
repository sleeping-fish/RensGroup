# @Author:Xinyu_Lu
# -*- coding = 'utf-8' -*-
# @Time:2020/9/18 15:54
# @File: transferTxt2Py.py
# @Software: PyCharm

from mpl_toolkits import mplot3d

import matplotlib.pyplot as plt
import numpy as np
x = np.load('x.npy')
y = np.load('y.npy')
z = np.load('z.npy').reshape(256, 256)

# 三维散点
# ax = plt.axes(projection='3d')
# ax.plot3D(x, y, z  ,)


r = np.linspace(0, 6, 20)
theta = np.linspace(-0.9 * np.pi, 0.8 * np.pi, 40)
r, theta = np.meshgrid(r, theta)
ax = plt.axes(projection='3d')
ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
plt.show()
