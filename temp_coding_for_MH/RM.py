import numpy as np
import math


# v1 = np.mat('-0.315708,0,0;0,-5.255654,0;0,0,-0.997695')
# v2 = np.mat('0.04019,0,0;0,1.619376,0;0,0,0.55102')
# v3 = np.mat('0.262933,0,0;0,-6.423854,0;0,0,1.41093')
# # print(v1, v2, v3)
# z = 18
# F = 1.469 * math.exp(-0.023 * math.pow(z,2)) + 2.679 * math.exp(-0.003 * math.pow((z-43.825),2))
#
# print(F)
# x = [0,  0.1 * F, F]
# E1 = np.mat(x)
# E2 = E1.T
# # print(E1 , E2)
#
# I1 = (E1 @ v1 @ E2) * (E1 @ v1 @ E2)
# I2 = (E1 @ v2 @ E2) * (E1 @ v2 @ E2)
# I3 = (E1 @ v3 @ E2) * (E1 @ v3 @ E2)
#
# print(I1, I2, I3)

def Angle(x, y):
    z = cp = np.cross(x, y)
    # print(z)
    # m = np.array(-1.42000008,  0.02999878,  0.33000183)
    # 分别计算两个向量的模：
    l_x = np.linalg.norm(x)
    l_y = np.linalg.norm(y)
    # print('向量的模=', l_x, l_y)

    # 计算两个向量的点积
    dian = x.dot(y)
    # print('向量的点积=', dian)

    # 计算夹角的cos值：
    cos_ = dian / (l_x * l_y)
    # print('夹角的cos值=', cos_)

    # 求得夹角（弧度制）：
    angle_hu = np.arccos(cos_)
    # print('夹角（弧度制）=', angle_hu)

    # 转换为角度值：
    angle_d = angle_hu * 180 / np.pi
    # print('夹角=%f°' % angle_d)

    return angle_hu


def RotateMat(x, y):  # 计算的是从y绕z轴转动到x的转动矩阵
    z = np.cross(x, y)
    z = z / np.linalg.norm(z)
    theta = Angle(x, y)
    c = math.cos(theta)
    s = math.sin(theta)
    # print(theta,z)
    temp = np.zeros((3, 3))  # 根据罗德里格方程计算转动矩阵
    temp[0][0] = c + (1 - c) * z[0] * z[0]
    temp[0][1] = (1 - c) * z[0] * z[1] - s * z[2]
    temp[0][2] = (1 - c) * z[0] * z[2] + s * z[1]
    temp[1][0] = (1 - c) * z[0] * z[1] + s * z[2]
    temp[1][1] = c + (1 - c) * z[1] * z[1]
    temp[1][2] = (1 - c) * z[1] * z[2] - s * z[0]
    temp[2][0] = (1 - c) * z[0] * z[2] - s * z[1]
    temp[2][1] = (1 - c) * z[1] * z[2] + s * z[0]
    temp[2][2] = c + (1 - c) * z[2] * z[2]
    # print(temp, c, s)
    return temp


def Rotate(x, y, z):  # 计算的是从y绕z轴转动到x的转动矩阵
    z = z / np.linalg.norm(z)
    x_v = (x @ z) * z  # 计算在转动轴上的投影1
    y_v = (x @ z) * z  # 计算在转动轴上的投影2
    m = x - x_v  # 计算垂直与转动轴的向量1
    n = y - y_v  # 计算垂直与转动轴的向量2
    theta = Angle(m, n)  # 计算二者夹角
    c = math.cos(theta)
    s = math.sin(theta)
    # print(theta, z, m, n)
    temp1 = np.zeros((3, 3))
    temp1[0][0] = c + (1 - c) * z[0] * z[0]
    temp1[0][1] = (1 - c) * z[0] * z[1] - s * z[2]
    temp1[0][2] = (1 - c) * z[0] * z[2] + s * z[1]
    temp1[1][0] = (1 - c) * z[0] * z[1] + s * z[2]
    temp1[1][1] = c + (1 - c) * z[1] * z[1]
    temp1[1][2] = (1 - c) * z[1] * z[2] - s * z[0]
    temp1[2][0] = (1 - c) * z[0] * z[2] - s * z[1]
    temp1[2][1] = (1 - c) * z[1] * z[2] + s * z[0]
    temp1[2][2] = c + (1 - c) * z[2] * z[2]
    # print(temp1, c, s)
    return temp1


def calculate(m, n, x, y):
    # m,n,x,y分别为：标准坐标a(主轴),读取坐标a'(主轴),标准坐标b,读取坐标b'
    # y: std ; x :transferred
    R = RotateMat(n, m)
    R1 = Rotate(y, x @ R, n)
    T = R @ R1
    return T


def SymMatrix(data):
    matrixSet = np.zeros((48, 3, 3))
    temp = np.zeros((3, 3))
    for i in range(data.shape[0]):
        temp[0][1] = data[i][1]
        temp[0][2] = data[i][2]
        temp[1][2] = data[i][4]
        tempDig = np.diag([data[i][0], data[i][3], data[i][5]])
        temp += tempDig
        temp += temp.T - np.diag(temp.diagonal())
        matrixSet[i] = temp
    return matrixSet


def test():
    m = np.array((2, 1, 0))  # 标准坐标a，主轴！
    n = np.array((1, 0, 2))  # 读取坐标a'，主轴！
    x = np.array((1, 1, 0))  # 标准坐标b
    y = np.array((1, 0, 1))  # 读取坐标b'
    T = calculate(m, n, x, y)
    print(T)


if __name__ == '__main__':
    test()
