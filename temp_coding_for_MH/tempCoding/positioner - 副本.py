"""
This programm reads in .pdb files, create a structure object to extract information about the position of certain atoms

More about the structure class, please refer to:
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc159
"""

import numpy as np
import pandas as pd
import re
import math
import time
from Bio.PDB.PDBParser import PDBParser


class RM:
    def __init__(self):
        pass

    @staticmethod
    def Angle(x, y):
        # z = cp = np.cross(x, y)
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

    @staticmethod
    def RotateMat(x, y):  # 计算的是从y绕z轴转动到x的转动矩阵
        z = np.cross(x, y)
        z = z / np.linalg.norm(z)
        theta = RM.Angle(x, y)
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

    @staticmethod
    def Rotate(x, y, z):  # 计算的是从y绕z轴转动到x的转动矩阵
        z = z / np.linalg.norm(z)
        x_v = (x @ z) * z  # 计算在转动轴上的投影1
        y_v = (x @ z) * z  # 计算在转动轴上的投影2
        m = x - x_v  # 计算垂直与转动轴的向量1
        n = y - y_v  # 计算垂直与转动轴的向量2
        theta = RM.Angle(m, n)  # 计算二者夹角
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

    @staticmethod
    def Calculate(m, n, x, y):
        # m,n,x,y分别为：标准坐标a(主轴),读取坐标a'(主轴),标准坐标b,读取坐标b'
        # y: std ; x :transferred
        R = RM.RotateMat(n, m)
        R1 = RM.Rotate(y, x @ R, n)
        T = R @ R1
        return T

    @staticmethod
    def SymMatrix(data):
        matrixSet = np.zeros((data.shape[0], 3, 3))
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


class Matrix:
    @staticmethod
    def LoadMartix(FilePath):
        df = pd.read_excel(FilePath)
        WaveNumber = df.WaveNumber.to_numpy()
        Label = df.label.to_numpy()
        Width = df.width.to_numpy()
        MatrixSet = RM.SymMatrix(df.to_numpy())
        return MatrixSet, Label, WaveNumber, Width

    @staticmethod
    def ProcessMartix(MatrixSet, WaveNumber, Label, Z, T_Set):
        S_Set = []
        WaveSet = []
        # 这边用字典存了T_set，不够可以加
        label_T_dict = {'1': T_Set[0], '2': T_Set[1], '3': T_Set[2], '4': T_Set[3]}
        for i in range(Label.shape[0]):
            T = label_T_dict[str(Label[i])]
            alpha = np.linalg.inv(T) @ MatrixSet[i] @ T
            s = np.array([0, 0, Z ** (-3)]).reshape(1, 3) \
                @ alpha \
                @ np.array([0, 0, Z ** (-3)]).reshape(3, 1)
            S_Set.append(s)
            WaveSet.append(WaveNumber[i])
        return S_Set, WaveSet

    @staticmethod
    def SaveData(DataSet, WaveNumber, Width):
        DataSet = (np.array(DataSet).reshape(-1, 1))
        WaveNumber = np.array(WaveNumber).reshape(-1, 1)
        Width = np.array(Width).reshape(-1, 1)
        temp = np.concatenate((WaveNumber, DataSet, Width), axis=1)
        return temp


# ------------------------------------------


# def read_single_Z(x):
#     x_Z = x[2]
#     return x_Z

# def read_Z(x, y):
#     x_Z = x[2]
#     y_Z = y[2]
#     print(abs((x_Z + y_Z) / 2))  # 可以以后加三次方用于估计场强作用


def Centre_Two_vector(x, y):
    m = x + y
    print(m[0:3])  # get_vector输出的第一位是中文，不是坐标矩阵


# def Tansform_mat(a,b,c,d):  #用于得到变化矩阵
#     m = np.mat(a)   # 标准朝向下的某一个向量矩阵a
#     n = np.mat(b)   # 标准朝向下的某一个向量矩阵b
#     o = np.mat(c)   # 实际坐标下的a对应的向量矩阵c
#     p = np.mat(d)   # 实际坐标下的b对应的向量矩阵d
#     e = np.dot(m.I, o) # 变化矩阵
#     l = np.dot(n, e)
#     if np.dot(p,l) / (np.linalg.norm(p) * np.linalg.norm(l)) == 1:  #判断夹角是否为0，为0则在两个向量均变化完成；如果不同，则需进一步变换（未完成，未验证）
#         return e
#     else:
#         e = np.dot(l.I, p)
#         return e
def SaveFiles(DataSet, SavePath):
    np.savetxt(SavePath, DataSet, fmt='%.20f')
    with open(SavePath, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write("@Author:Xinyu_Lu\n@Time:{}\n@Total Lines:{}\n"
                .format(time.asctime(time.localtime(time.time())),
                        DataSet.shape[0]) + content)
        f.close()


class Protein:

    def __init__(self, filename, MatrixDir, TxtSaveDir, structure_id):
        self.filename = filename
        self.structure_id = structure_id
        self.MatrixDir = MatrixDir
        self.TxtSaveDir = TxtSaveDir

    def toStructure(self):
        parser = PDBParser(PERMISSIVE=1)
        self.structure = parser.get_structure(self.structure_id, self.filename)

    def getResidue(self, SaveName):
        # f = open("POI.txt", "+a")

        # structure contains model?
        if len(self.structure) != 1:
            return
        else:
            self.model = self.structure[0]
        # model contains chain?
        if len(self.model) != 1:
            return
        else:
            self.chain = self.model["A"]
        # print residues (sequence)
        self.aaList = self.chain.get_list()
        print('组成该蛋白的氨基酸个数：', len(self.aaList),
              '\nfile=%s' % self.filename)
        # print('Amino acid:','        vector1:', '     vector2:', '    vector3:', '      Distance from substrate',
        # 'file=f')
        for aa in self.aaList:
            atoms = list(aa.get_atoms())  # 将氨基酸残基内的原子收集
            # print(aa.get_resname())
            Data_TRP = np.zeros((1, 3))
            Data_HIS = np.zeros((1, 3))
            Data_ASP = np.zeros((1, 3))
            Data_ALA = np.zeros((1, 3))
            Data_GLU = np.zeros((1, 3))
            Data_PHE = np.zeros((1, 3))
            Data_ARG = np.zeros((1, 3))
            Data_SER = np.zeros((1, 3))
            Data_GLY = np.zeros((1, 3))
            Data_TYR = np.zeros((1, 3))
            Data_VAL = np.zeros((1, 3))
            if aa.get_resname() == "HIS":  # 如果残基为组氨酸，则采一部分原子，需要5个向量
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[6].get_vector()  # atom：CD2
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CB
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[8].get_vector()  # atom: C
                vector_1 = Atom6 - Atom1  # 酰胺键轴方向c-c(o), [-0.751, 1.31, -0.19]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.97, -1.083, -0,131]
                vector_3 = Atom4 - Atom1  # c-c,用于ch2, [1.175, -0.224, -0.996]
                vector_4 = Atom4 - Atom3  # 咪唑环向量1, [-1.344, -0.488, -0.436]
                vector_5 = Atom3 - Atom2  # 咪唑环向量2, [-1.159, 0.684, -0.284]
                cp = np.cross(vector_4[0:3], vector_5[0:3])  # 咪唑环平面法向量
                Centre_Two_vector(vector_4, vector_1)  # 得到某两个向量的角平分线，主要用于CH2的振动计算

                # print(vector_5)
                # print(vector_5[1:3])
                # print(aa.get_resname(), cp, vector_4[0:3], vector_5[0:3], '    CG原子的z轴坐标:', Atom3[2], 'file=f')
                # read_Z(Atom1, Atom3)  #通过计算两个原子的z坐标，估计分子距离基底的高度

            if aa.get_resname() == "TYR":  # 如果残基为酪氨酸，则采一部分原子
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[10].get_vector()  # tyr_atom: C
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [0.066, -1.528, 0.013]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [1.306, 0.59, -0.309]
                vector_3 = Atom2 - Atom1  # for the future, [-1.037, 0.484, -1.033]
                vector_4 = Atom3 - Atom2  # 苯环向量1, [-1.418, -0.44, 0.288]
                vector_5 = Atom3 - Atom4  # 苯环向量2, [0.556, 1.123, 0.619]
                # cp = np.cross(vector_4[0:3], vector_5[0:3])  # 咪唑环平面法向量
                # print(aa.get_resname(), cp, vector_4[0:3], vector_5[0:3],'    CG原子的z轴坐标:', Atom3[2], 'file=f')

            if aa.get_resname() == "PHE":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[9].get_vector()  # atom: C
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-0.487, 1.424, 0.233]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-1.137, -0.799, -0.45]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.151, -0.118, -1.036]
                vector_4 = Atom3 - Atom2  # 苯环向量1, [1.384, 0.386, 0.474]
                vector_5 = Atom3 - Atom4  # 苯环向量2, [0.686, 1.078, -0.567]
                # cp = np.cross(vector_4[0:3], vector_5[0:3])  # 咪唑环平面法向量
                # print(aa.get_resname(), cp, vector_4[0:3], vector_5[0:3], '    CG原子的z轴坐标:', Atom3[2], 'file=f')

            if aa.get_resname() == "TRP":
                '''计算Atoms 和 vector'''
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[11].get_vector()  # atom: CD2
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[12].get_vector()  # atom: C
                vector_1 = (Atom6 - Atom1)[0: 3]  # 键轴方向c-c(o), [-0.515, 1.208, 0.771]
                vector_2 = (Atom5 - Atom1)[0: 3]  # c-n,暂时没用, [-1.18, -0.781, -0.38]
                vector_3 = (Atom2 - Atom1)[0: 3]  # for the future,暂时没用,[0.88, 0.472, -1.189]
                vector_4 = (Atom4 - Atom3)[0: 3]  # 苯环向量1, [1.23, -0.564, 0.501]
                vector_5 = (Atom3 - Atom2)[0: 3]  # 苯环向量2, [1.241, 0.745, 0.395]
                # cp = np.cross(vector_4[0:3], vector_5[0:3])  # 咪唑环平面法向量

                '''计算不同的T '''
                # 苯环标准向量1 苯环标准向量2
                phStd1 = np.array([1.23, -0.564, 0.501])
                phStd2 = np.array([1.241, 0.745, 0.395])
                # 苯环T记为1
                T1 = RM.Calculate(phStd1, vector_4, phStd2, vector_5)
                # c - c(o)
                CCO_Std = np.array([-0.515, 1.208, 0.771])
                # c-n
                CN_Std = np.array([-1.18, -0.781, -0.38])
                T2 = RM.Calculate(CCO_Std, vector_1, CN_Std, vector_2)
                # 这边T3和T4可以再写，暂时用T1和T2占位
                T3 = T1
                T4 = T2
                T_Set = [T1, T2, T3, T4]
                Z = Atom3[2]

                '''处理Matrix'''
                MatrixSet, Label, WaveNumber, Width = Matrix.LoadMartix(
                    FilePath=self.MatrixDir + fr'\{aa.get_resname()}.xlsx')
                S_Set, WaveSet = Matrix.ProcessMartix(MatrixSet, WaveNumber, Label, Z, T_Set)
                temp = Matrix.SaveData(S_Set, WaveSet, Width)
                np.append(Data_TRP, temp)
            if aa.get_resname() == "THR":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG2
                Atom4 = atoms[4].get_vector()  # atom: OG1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[5].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: O
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-1.38, -0.451, 0.466]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.156, 1.407, 0.32]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.067, -0.976, 0.568]
                vector_4 = Atom3 - Atom2  # c-c, [1.406, 0.292, -0.517]
                vector_5 = Atom4 - Atom2  # c-o, [-0.076, 0.116, 1.423]

            if aa.get_resname() == "ALA":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                # Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: C
                Atom5 = atoms[0].get_vector()  # atom: N
                # Atom6 = atoms[10].get_vector()  # atom: C
                vector_1 = Atom4 - Atom1  # 键轴方向c-c(o), [1.474, -0.176, 0.33]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.719, -1.194, 0.434]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-0.611, 1.255, 0.634]
                # vector_4 = Atom3 - Atom2  # 苯环向量1
                # vector_5 = Atom3 - Atom4  # 苯环向量2

            if aa.get_resname() == "ILE":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG2
                Atom4 = atoms[4].get_vector()  # atom: CG1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # C
                Atom7 = atoms[5].get_vector()  # atom: CD
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [0.898, 1.084, 0.581]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.674, -1.288, 0.162]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-1.424, -0.07, 0.635]
                vector_4 = Atom3 - Atom2  # c-c 1, [0.072, -0.046, 1.53]
                vector_5 = Atom4 - Atom2  # c-c 2, [-0.964, 1.071, -0.539]
                vector_6 = Atom7 - Atom4  # c-c 3, [-1.466, -0.242, 0.374]

            if aa.get_resname() == "LEU":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # C
                Atom7 = atoms[5].get_vector()  # atom: CD2
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [0.937, -1.201, -0.028]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.728, 1.143, -0.551]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-1.319, -0.207, -0.78]
                vector_4 = Atom4 - Atom3  # c-c 1, [-0.575, 0.838, 1.152]
                vector_5 = Atom7 - Atom3  # c-c 2, [-1.107, -0.362, -0.999]

            if aa.get_resname() == "VAL":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG1
                Atom4 = atoms[4].get_vector()  # atom: CG2
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[5].get_vector()  # C
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-1.004, -0.966, 0.614]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.484, 1.363, 0.212]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.45, -0.135, 0.559]
                vector_4 = Atom3 - Atom2  # c-c 1, [0.03, 0.032, 1.532]
                vector_5 = Atom4 - Atom2  # c-c 2, [-.778, -1.19, -0.571]

            if aa.get_resname() == "LYS":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[7].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: CE
                Atom8 = atoms[6].get_vector()  # atom: NZ
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-0.709, 1.304, -0.348]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.877, -1.1, -0.395]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.386, -0.174, -0.659]
                vector_4 = Atom3 - Atom2  # c-c 1, [1.161, 0.719, 0.693]
                vector_5 = Atom3 - Atom4  # c-c 2, [-1.345, 0.283, 0.685]
                vector_6 = Atom7 - Atom4  # c-c 3, [1.135, 0.715, 0.728]
                vector_7 = Atom7 - Atom8  # c-n 2, [-1.348, 0.236, 0.637]

            ###########################

            if aa.get_resname() == "ARG":
                Atom1 = atoms[1].get_vector()  # atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[9].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: NE
                Atom8 = atoms[6].get_vector()  # atom: CZ
                Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [0.75, 1.296, 0.281]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.874, -1.107, 0.281]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-1.36, -0.131, 0.72]
                vector_4 = Atom3 - Atom2  # c-c 1, [-1.178, 0,708, 0.043]
                vector_5 = Atom3 - Atom4  # c-c 2, [1.314, 0.239, -0.742]
                vector_6 = Atom7 - Atom4  # c-n 1, [-1.115, 0.684, -0.742]
                vector_7 = Atom7 - Atom8  # c-n 2, [1.318, 0.31, -0.404]
                vector_8 = Atom9 - Atom8  # c-n 3, [-0.358, -1.334, -0.269]

            if aa.get_resname() == "ASN":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: OD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: ND2
                # Atom8 = atoms[6].get_vector()  # atom: NZ
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [0.668, -1.363, -0.087]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.951, 1.007, -0.472]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-1.316, 0.094, -0.791]
                vector_4 = Atom3 - Atom2  # c-c 1, [-1.193, -0.595, 0.743]
                vector_5 = Atom3 - Atom4  # c-o 2, [0.167, -0.173, -1.195]
                vector_6 = Atom3 - Atom7  # c-n 3, [0.908, 0.66, 0.776]
                # vector_7 = Atom7 - Atom8  # c-c 3

            if aa.get_resname() == "ASP":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: OD1
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: OD2
                # Atom8 = atoms[6].get_vector()  # atom: NZ
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [0.947, 1.141, 0.389]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.651, -1.297, -0.017]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-1.238,-0.047 ,0.948]
                vector_4 = Atom3 - Atom2  # c-c 1, [-1.275, 0.698, -0.482]
                vector_5 = Atom3 - Atom4  # c-o 2, [-0.112, -1.284, 0.403]
                vector_6 = Atom3 - Atom7  # c-o 3, [1.066, 0.552, 0.008]
                # vector_7 = Atom7 - Atom8  # c-c 3

            if aa.get_resname() == "CYS":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: s
                Atom4 = atoms[4].get_vector()  # atom: C
                Atom5 = atoms[0].get_vector()  # atom: N
                # Atom6 = atoms[4].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: CE
                # Atom8 = atoms[6].get_vector()  # atom: NZ
                vector_1 = Atom4 - Atom1  # 键轴方向c-c(o), [-0.691, -1.243, 0.556]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.96, 1.096, 0.133]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.329, 0.245, 0.737]
                vector_4 = Atom3 - Atom2  # c-s 1, [1.224, -1.337, -0.375]
                # vector_5 = Atom3 - Atom4  # c-c 2
                # vector_6 = Atom7 - Atom4  # c-c 3
                # vector_7 = Atom7 - Atom8  # c-c 3

            if aa.get_resname() == "GLN":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[7].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: OE1
                Atom8 = atoms[6].get_vector()  # atom: NE2
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-0.842, 1.215, -0.348]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.706, -1.187, -0.459]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.43, 0.002, -0.599]
                vector_4 = Atom3 - Atom2  # c-c 1, [1.013, 0.914, 0.695]
                vector_5 = Atom3 - Atom4  # c-c 2, [-1.432, 0.335, 0.412]
                vector_6 = Atom7 - Atom4  # c-o 1, [0.4, -1.146, -0.117]
                vector_7 = Atom8 - Atom4  # c-n 2, [0.727, 1.076, -0.2]
                # vector_8 = Atom9 - Atom8  # c-n 3

            if aa.get_resname() == "GLU":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[7].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: OE1
                Atom8 = atoms[6].get_vector()  # atom: OE2
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-0.84, -1.185, 0.466]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.71, 1.221, 0.357]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.425, 0.049, 0.599]
                vector_4 = Atom3 - Atom2  # c-c 1, [0.999, -0.988, -0.613]
                vector_5 = Atom3 - Atom4  # c-c 2, [-1.413, -0.253, -0.613]
                vector_6 = Atom7 - Atom4  # c-o 1, [0.923, -0.854, -0.105]
                vector_7 = Atom8 - Atom4  # c-o 2, [0.425, 1.078, 0.484]
                # vector_8 = Atom9 - Atom8  # c-n 3

            if aa.get_resname() == "GLY":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                # Atom2 = atoms[2].get_vector()  # atom：CB
                # Atom3 = atoms[3].get_vector()  # atom: CG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[2].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: NE
                # Atom8 = atoms[6].get_vector()  # atom: CZ
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [1.268, 0.814, -0.114]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-1.157, 0.811, -0.342]
                # vector_3 = Atom2 - Atom1  # for the future,暂时没用
                # vector_4 = Atom3 - Atom2  # c-c 1
                # vector_5 = Atom3 - Atom4  # c-c 2
                # vector_6 = Atom7 - Atom4  # c-n 1
                # vector_7 = Atom7 - Atom8  # c-n 2
                # vector_8 = Atom9 - Atom8  # c-n 3

            if aa.get_resname() == "MET":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: CG
                Atom4 = atoms[4].get_vector()  # atom: SD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[6].get_vector()  # atom: C
                Atom7 = atoms[5].get_vector()  # atom: CE
                # Atom8 = atoms[6].get_vector()  # atom: CZ
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [-0.738, 1.263, -0.426]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [-0.857, -1.141, -0.308]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [1.379, -0.195, -0.672]
                vector_4 = Atom3 - Atom2  # c-c 1, [1.142, 0.81, 0.613]
                vector_5 = Atom3 - Atom4  # c-s 1, [-1.583, 0.456, 0.812]
                vector_6 = Atom7 - Atom4  # c-s 2, [1.123, 1.133, 0.884]
                # vector_7 = Atom7 - Atom8  # c-n 2
                # vector_8 = Atom9 - Atom8  # c-n 3

            if aa.get_resname() == "SER":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CA
                Atom2 = atoms[2].get_vector()  # atom：CB
                Atom3 = atoms[3].get_vector()  # atom: OG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[4].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: NE
                # Atom8 = atoms[6].get_vector()  # atom: CZ
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom1  # 键轴方向c-c(o), [1.187, -0.789, -0.556]
                vector_2 = Atom5 - Atom1  # c-n,暂时没用, [0.116, 1.434, -0.263]
                vector_3 = Atom2 - Atom1  # for the future,暂时没用, [-1.293, -0.513, -0.631]
                vector_4 = Atom3 - Atom2  # c-c 1, [-0.105, -1.39, 0.294]
                # vector_5 = Atom3 - Atom4  # c-c 2
                # vector_6 = Atom7 - Atom4  # c-n 1
                # vector_7 = Atom7 - Atom8  # c-n 2
                # vector_8 = Atom9 - Atom8  # c-n 3

            if aa.get_resname() == "PRO":
                Atom1 = atoms[1].get_vector()
                # print(atoms[1])  #  atom：CD
                Atom2 = atoms[4].get_vector()  # atom：CA
                # Atom3 = atoms[3].get_vector()  # atom: CG
                # Atom4 = atoms[4].get_vector()  # atom: CD
                Atom5 = atoms[0].get_vector()  # atom: N
                Atom6 = atoms[5].get_vector()  # atom: C
                # Atom7 = atoms[5].get_vector()  # atom: NE
                # Atom8 = atoms[6].get_vector()  # atom: CZ
                # Atom9 = atoms[7].get_vector()  # atom: NH1
                vector_1 = Atom6 - Atom2  # 键轴方向c-c(o), [1.352, 0.315, -0.613]
                vector_2 = Atom5 - Atom1  # c-n, 环1, [1.136, 0.21, 0.912]
                vector_3 = Atom2 - Atom1  # 环2, [2.024, -0.955, 0.853]
                # vector_4 = Atom3 - Atom2  # c-c 1
                # vector_5 = Atom3 - Atom4  # c-c 2
                # vector_6 = Atom7 - Atom4  # c-n 1
                # vector_7 = Atom7 - Atom8  # c-n 2
                # vector_8 = Atom9 - Atom8  # c-n 3

            '''在这个地方加几行，把上面每个if返回的data合并，然后写入文件'''
            '''我的想法是，在每个if下加一个SaveData函数，保存计算出来的data，最后再拼接起来'''

            DataSet = np.concatenate((Data_TRP,
                                      Data_HIS,
                                      Data_ASP,
                                      Data_ALA,
                                      Data_GLU,
                                      Data_PHE,
                                      Data_ARG,
                                      Data_SER,
                                      Data_GLY,
                                      Data_TYR,
                                      Data_VAL), axis=0)
            DataSet = DataSet[np.all(DataSet != 0, axis=1)]
            SaveFiles(DataSet=DataSet, SavePath=self.TxtSaveDir + '\\' + SaveName)

    # def getPosition(self):
    #     for aa in self.aaList:
    #         if aa == "PHE" or aa == "TYR":
    #             print(aa.get_resname)
    #             atomList = aa.get_atoms()
    #             positionAtom1 = atomList[0].get_pos()
    #             positionAtom2 = atomList[0].get_pos()
    #             positionAtom3 = atomList[0].get_pos()
    #             getNormVector(positionAtom1, positionAtom2, positionAtom3)
    #         elif aa.get_resname() == "TRP":
    #             atomList = aa.get_atoms()
    #             positionAtom1 = atomList[0].get_pos()
    #             positionAtom2 = atomList[0].get_pos()
    #             positionAtom3 = atomList[0].get_pos()
    #             getNormVector(positionAtom1, positionAtom2, positionAtom3)
    #
    # def getNormVector(a, b, c):
    #     vector1 = b - a
    #     vector2 = b - c
    #     cp = np.cross(vector1, vector2)
