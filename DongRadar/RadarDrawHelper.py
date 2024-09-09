# -*- coding:utf-8 -*-
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# plt.switch_backend('agg')

def radar_colormap(colortable_type="Z", proportion=True, sep=False, spacing='c'):
    '''
    根据给定的numpy_array矩阵，返回colormap。

    numpy_array: string
        numpy_array矩阵,结构如下: [[value, R, G, B],[]]
    proportion: boolean, default is True 比例，使得均衡
        proportion为True时LinearSegmentedColormap将由输入值给出
    sep: boolean, default is False
        当sep为真时，colormap将设置为颜色块，没有颜色梯度变化
    spacing: string, 默认为'c', 仅仅当sep为False时可用
        When spacing is 'c', the color blocks will be equally spaced.
        A ListedColormap will be returned.颜色块将被等距隔开
        When spacing is 'v', the length of color blocks will be based
        on the input values. A LinearSegmentedColormap will be returned.
        颜色块的长度将基于输入值
    '''
    import numpy as np
    import matplotlib.colors as cmx
    if colortable_type not in ["dBZ","W","V"]:
        raise ValueError("colortable_type not in [dBZ,W,V]")
    #中国局雷达数据颜色表
    arr = None
    if colortable_type=="dBZ":
        arr = np.array([
            [0,255,255,255],
            [5,192,192,254],
            [10,122,114,238],
            [15,30, 38, 208],
            [20,166, 252,168],
            [25,0, 234, 0],
            [30,16, 146, 26],
            [35,252,244,100],
            [40,200, 200, 2],
            [45,140,140,0],
            [50,254,172,172],
            [55,254,100,84],
            [60,238,2,48],
            [65,212,142,254],
            [70,170,36,250]
            ])
    if colortable_type=="W":
        arr = np.array([
            [0,255,255,255],
            [0.04,0,107,253],
            [0.08,0,186,253],
            [0.12,111,248,255],
            [0.16,0,150,50],
            [0.20,0,220,0],
            [0.23,180,255,180],
            [0.26,196,166,0],
            [0.29,255,255,0],
            [0.32,238,255,0],
            [0.35,255,0,0],
            [0.38,255,100,100],
            [0.41,255,180,180],
            [0.44,150,0,180],
            [0.47,200,100,155],
            [0.50,341,98,153]
            ])
    if colortable_type=="V":
        arr = np.array([
            [-30,126,224,254],
            [-27,0,224,254],
            [-20,0,176,176],
            [-15,0,254,0],
            [-10,0,196,0],
            [-5,0,128,0],
            [-1,254,254,254],
            [0,252,252,252],
            [1,254,0,0],
            [5,254,88,88],
            [10,254,176,176],
            [15,254,124,0],
            [20,255,210,0],
            [27,254,254,0]
            ])
    inidict = {'red':None, 'green':None, 'blue':None}
    if sep == True:
        if spacing == 'c':
            ar1 = arr.transpose()
            value = []
            count = 0
            while count < len(ar1[1]):
                value.append((int(ar1[1][count]) / 255, int(ar1[2][count]) / 255,
                              int(ar1[3][count]) / 255))
                count = count + 1
            return cmx.ListedColormap(value, 256)
        elif spacing == 'v':
            value = []
            r = []
            g = []
            b = []
            for i in arr:
                value.append(float(i[0]))
                r.append(int(i[1]) / 255)
                g.append(int(i[2]) / 255)
                b.append(int(i[3]) / 255)
                if len(value) > 1:
                    if value[-1] < value[-2]:
                        raise ValueError('Values must be in order')
            inivalue = value[0]
            maxvalue = value[-1]
            drange = maxvalue - inivalue
            rpart = []
            gpart = []
            bpart = []
            count = 0
            try:
                while count < len(value):
                    tupr = ((value[count] - inivalue) / drange, r[count], r[count])
                    tupg = ((value[count] - inivalue) / drange, g[count], g[count])
                    tupb = ((value[count] - inivalue) / drange, b[count], b[count])
                    tupra = ((value[count + 1] - inivalue-1e-5) / drange, r[count], r[count])
                    tupga = ((value[count + 1] - inivalue-1e-5) / drange, g[count], g[count])
                    tupba = ((value[count + 1] - inivalue-1e-5) / drange, b[count], b[count])
                    rpart.append(tupr)
                    rpart.append(tupra)
                    gpart.append(tupg)
                    gpart.append(tupga)
                    bpart.append(tupb)
                    bpart.append(tupba)
                    count=count+1
            except IndexError:
                rpart.append((1, r[count], r[count]))
                gpart.append((1, g[count], g[count]))
                bpart.append((1, b[count], b[count]))
            inidict['red'] = rpart
            inidict['green'] = gpart
            inidict['blue'] = bpart
            return cmx.LinearSegmentedColormap('my_colormap', inidict, 256)

    elif sep == False:
        value = []
        r = []
        g = []
        b = []
        for i in arr:
            value.append(float(i[0]))
            r.append(int(i[1]) / 255)
            g.append(int(i[2]) / 255)
            b.append(int(i[3]) / 255)
            if len(value) > 1:
                if value[-1] < value[-2]:
                    raise ValueError('Values must be in order')
        inivalue = value[0]
        maxvalue = value[-1]
        drange = maxvalue-inivalue
        rpart = []
        gpart = []
        bpart = []
        if proportion == True:
            count = 0
            while count < len(value):
                tupr = ((value[count] - inivalue) / drange, r[count], r[count])
                tupg = ((value[count] - inivalue) / drange, g[count], g[count])
                tupb = ((value[count] - inivalue) / drange, b[count], b[count])
                rpart.append(tupr)
                gpart.append(tupg)
                bpart.append(tupb)
                count=count + 1
        elif proportion == False:
            count = 0
            while count < len(value):
                tupr = (count / len(value), r[count], r[count])
                tupg = (count / len(value), g[count], g[count])
                tupb = (count / len(value), b[count], b[count])
                rpart.append(tupr)
                gpart.append(tupg)
                bpart.append(tupb)
                count = count + 1
        inidict['red'] = rpart
        inidict['green'] = gpart
        inidict['blue'] = bpart
        return cmx.LinearSegmentedColormap('my_colormap', inidict, 256)