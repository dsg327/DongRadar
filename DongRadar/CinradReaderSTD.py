# -*- coding: utf-8 -*-
# ===========================================================================
# Filename: CinradReaderSTD.py
# Description: Read of Radar Basic Data Files
#
# Author: Dong Hongchang
# Email: dsg327@163.com
# License: MIT License
#
# Copyright (c) 2023 DongHc
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ===========================================================================

"""
Read of Radar Basic Data Files

This script can be used to:
- Read of Radar Basic Data Files
- Provide radar data query services

Usage:
    reader=CinradReaderSTD(filename)
"""

import numpy as np
import struct
from pathlib import Path


class RadarError(Exception):
    def __init__(self, description):
        self.dsc = description
    def __str__(self):
        return repr(self.dsc)
    

class CinradReaderSTD():
    '''标准格式雷达数据读取类'''

    def __init__(self, filepath=None):
        self.last_error=''
        self.data = {
            'common_block':{
                'generic_header':{},
                'site_conf':{},
                'task_conf':{},
                'cut_conf':[]
            }, 
            'radial_block':[],
            "radial": {}
        }
        self.vcp_standard_elevation={
            'VCP11':[0.50,1.45,2.40,3.35,4.30,5.25,6.2,7.5,8.7,10.00,12.00,14.00,16.70,19.50],
            'VCP11D':[0.48,0.48,1.49,1.49,2.42,3.34,4.31,5.23,6.20,7.51,8.70,10.02,12.00,14.02,16.70,19.51],
            'VCP21':[0.50,1.45,2.40,3.35,4.30,6.00,9.90,14.6,19.5],
            'VCP21D':[0.53,0.53,1.45,1.45,2.42,3.43,4.31,6.02,9.93,14.59,19.51],
            'VCP31':[0.50,1.5,2.50,3.5,4.50],
            'VCP31D':[0.48,0.48,1.49,1.49,2.50,2.50,3.52,4.48],
            'VCP32':[0.50,1.5,2.50,3.5,4.50],
            'VCP32D':[0.48,0.48,1.49,1.49,2.50,3.52,4.48]
        }
        self.filepath = filepath
        if filepath!=None:
            self.read_data(filepath)
            
    def get_error(self):
        return self.last_error
            
    def get_data(self):
        return self.data
    
    def get_all_datatype(self):
        '''获取所有的数据类型'''
        return self.data['radial'].keys()
    
    def get_unique_radial_item(self,datatype='dBZ',item='elevation'):
        '''
        获取所有扫描射线的不重复数据项(属性)
        
        Parameters:
        datatype (str): 数据类型，例如dBZ、ZDR、V、W等
        item (str): 数据项目，可选的例如elevation、ppiElevation等
        '''
        return list(set(self.data['radial'][datatype][item]))
    
    def get_values_bylatlon(self, lat, lon, datatype="dBZ"):
        """
        通过经纬度获取所有层的回波数据值
        """
        result={'height':[],'value':[]}
        for layerid,v in enumerate(self.data['common_block']['cut_conf']):
            if not self.has_element(layerid,datatype):
                continue
            height,value = self.get_value_bylatlon(lat,lon,layerid,datatype)
            result['height'].append(height)
            result['value'].append(value)
        return result
    
    def get_value_bylatlon(self, lat, lon, layerid, datatype="dBZ"):
        """
        通过经纬度获取特定层的回波数据值
        """
        site_lat=self.data['common_block']['site_conf']['latitude']
        site_lon=self.data['common_block']['site_conf']['longitude']
        site_alt=self.data['common_block']['site_conf']['antennaHeight']

        from .RadarCalculationHelper import haversine_distance
        from .RadarCalculationHelper import calculate_bearing

        distance=haversine_distance(lat,lon,site_lat,site_lon)*1000
        azimuth=calculate_bearing(site_lat,site_lon,lat,lon)
        # print(lat,lon,site_lat,site_lon)
        # print(distance, azimuth)

        return self.get_value_bypolar(distance, azimuth, layerid, datatype)

    def get_value_bypolar(self, distance, azimuth, layerid, datatype="dBZ"):
        """
        通过距雷达的地面直线距离distance与方位角θ获取特定层的雷达数据
        
        Parameters:
        distance为目标点距离雷达站的直线距离
        azimuth 方位角
        雷达体扫层编号
        数据类型

        Return:
        height,value
        """
        # 获取数据
        start_pos,end_pos=self.get_cut_data_index(layerid,datatype)
        radar_data=self.data['radial'][datatype]['data'][start_pos:end_pos]
        azimuth_list=self.data['radial'][datatype]['azimuth'][start_pos:end_pos]

        layer_elevation=self.get_standard_elevation_byid(layerid,self.data['common_block']['task_conf']['taskName'])
        radius = distance / np.cos(np.deg2rad(layer_elevation))
        
        ku_count = len(radar_data[0])
        ku_length = self.get_ku_length(layerid,datatype)
        max_radius = ku_length * ku_count
        if radius > max_radius:
            # print("radius：%s超出了雷达的最大可探测范围%s！" % (radius, max_radius))
            return None
        
        diff_azimuth_arr = abs(np.array(azimuth_list) - azimuth)
        azimuth_index = np.where(diff_azimuth_arr == np.min(diff_azimuth_arr))[0][0]
        # print("查找到的方位角索引值为%s,对应的θ值为%s,ρ值为%s,库长为%s,数据索引值为%s."%
        #     (find_theta_index, find_theta,rho,ku_length,int(rho / ku_length)))
        # print("你查询点的数据值为：%s"%target_data[find_theta_index][int(rho / ku_length)])

        # 计算离地高度, 要考虑到地球曲率
        height = np.sin(np.deg2rad(layer_elevation))*distance
        r = 6371393 # 地球半径
        rm = 8500000 # 等效地球半径
        a=np.sqrt(2*(r**2)-2*(r**2)*np.cos(distance/r))
        b=np.sqrt(rm**2+a**2-2*rm*a*np.cos((np.pi-distance/r)/2))
        h=rm-b
        # h=r/np.cos(distance/r)-r
        # h=h*np.cos(distance/r)
        # print('height is %f,h is %f'%(height,h))
        return (height+h,radar_data[azimuth_index][int(radius / ku_length)])
        
    
    def get_ppi_data(self,layerid,show_range=330000,datatype='dBZ'):
        """
        获取PPI数据
        
        Parameters:
        layerid (float): 显示哪层,vcp layer,从1开始编号
        show_range (float): 最大显示范围，单位米
        datatype (str): PPI
        
        Return:
        (radius,azimuth,elevation,show_data)
        """
        start_pos,end_pos=self.get_cut_data_index(layerid,datatype)
        # 获取连续的数据
        radar_data=self.data['radial'][datatype]['data'][start_pos:end_pos]
        dataLength=self.data['radial'][datatype]['dataLength'][start_pos:end_pos]
        lentypelist=list(set(dataLength))
        if len(lentypelist)!=1:
            self.last_error='检测到一个ppi屏幕的不同径向数据库数不同。'
            return None
        radar_data=np.array(radar_data).astype(np.float64)
        azimuth=self.data['radial'][datatype]['azimuth'][start_pos:end_pos]
        elevation=self.data['radial'][datatype]['elevation'][start_pos:end_pos]
        # 裁剪雷达数据
        max_range = 0
        ku_length = 0
        ku_count = lentypelist[0]
        ku_length = self.get_ku_length(layerid,datatype)
        max_range = ku_length * ku_count
        if show_range > max_range:
            print("传入的show_range：%s超出了雷达的最大可探测范围%s，系统将自动进行重置！" % (show_range, max_range))
            show_range = max_range
        show_data=radar_data.transpose()[:int(show_range / ku_length)].transpose()
        # print(show_data.shape)
        show_data= np.ma.array (show_data, mask=np.isnan(show_data))#nan数据标记为bad
        # 获得极轴坐标轴
        # rho = np.arange(ku_length, int(show_range) + ku_length, ku_length) #该方法在显示范围不是库长整数倍时rho坐标轴比数据多一个
        radius = list(range(ku_length, ku_length+ku_length*show_data.shape[1], ku_length))
        return (radius,azimuth,elevation,show_data)
    
    def get_cut_data_index(self,cutid,datatype='dBZ'):
        """
        获取一次扫描切片的数据位置
        """
        # 获取连续的数据
        if not self.has_element(cutid,datatype):
            raise RadarError("当前要素在当前层不存在!")
        # 找到数据的起始和结束位置
        start_pos,end_pos = None,None
        continuous_state=False
        for i, value in enumerate(self.data['radial'][datatype]['elevationNumber']):
            if continuous_state==False:
                if value == cutid:
                    start_pos = i
                    continuous_state=True
            else:
                if value !=cutid:
                    end_pos = i
                    break
        return (start_pos,end_pos)
    
    def has_element(self,cutid,datatype='dBZ'):
        """
        返回cut中是否存在某个要素
        """
        if cutid in self.data['radial'][datatype]['elevationNumber']:
            return True
        else:
            return False
    
    def get_ku_length(self,layerid,datatype='dBZ'):
        """
        获取PPI扫描中某一层的库长
        """
        if datatype in ['V','W']:
            ku_length = self.data['common_block']['cut_conf'][layerid-1]['dopplerResolution']
        else:
            ku_length = self.data['common_block']['cut_conf'][layerid-1]['logResolution']
        return ku_length


        
    def get_datatype_byid(self,typeid):
        '''通过数据类型id获得气象数据类型key'''
        type_dic={
            1:'dBT',
            2:'dBZ',
            3:'V',
            4:'W',
            5:'SQI',
            6:'CPA',
            7:'ZDR',
            8:'LDR',
            9:'CC',
            10:'DP',
            11:'KDP',
            12:'CP',
            14:'HCL',
            15:'CF',
            16:'SNRH',
            17:'SNRV',
            32:'Zc',
            33:'Vc',
            34:'Wc',
            35:'ZDRc'
        }
        if typeid not in type_dic.keys():
            raise RadarError("数据类型代码%d不存在，请核实!"%typeid)
        return type_dic[typeid]
        
    def rhexgbk(self,bytes_data,encoding=None):
        """从二进制block_data读取slen长度的字符并用encoding解码，以x00作为字符串结束"""
        #str(struct.unpack("16s",block_data[0:16])[0],'gb2312').rstrip("\x00")有可能0x00后乱码字符
        if encoding==None:
            return bytes_data[:bytes_data.find(0x00)]
        zs=bytes_data.find(0x00)
        if zs==-1:
            return str(bytes_data,encoding)
        return str(bytes_data[:bytes_data.find(0x00)],encoding)

    def read_data(self, filepath=None):
        '''读取并解析文件'''
        path = Path(filepath)
        filename = path.name
        filetype = path.suffix
        # 开始读取文件数据
        f = None
        if filetype.endswith('bz2'):
            import bz2
            f = bz2.open(filepath, 'rb')
        else:
            f = open(filepath, 'rb')
        # 通用头
        block_data=f.read(32)
        pdata={}
        pdata["magicNumber"]=self.rhexgbk(struct.unpack("4s",block_data[0:4])[0],'gb2312') #固定标识(RSTM)
        if pdata["magicNumber"]!='RSTM':
            raise RadarError("雷达文件类型并非STD，固定标识检查不一致，请核实!")
        pdata["majorVersion"]=struct.unpack("H",block_data[4:6])[0] #主版本号
        pdata["minorVersion"]=struct.unpack("H",block_data[6:8])[0] #次版本号
        pdata["genericType"]=struct.unpack("i",block_data[8:12])[0] #文件类型
        pdata["productType"]=struct.unpack("i",block_data[12:16])[0] #产品类型
        pdata["reserved"]=self.rhexgbk(struct.unpack("16s",block_data[16:32])[0]) #保留字段
        self.data['common_block']['generic_header']=pdata
        # 站点配置信息
        block_data=f.read(128)
        pdata={}
        pdata['siteCode']=self.rhexgbk(struct.unpack("8s",block_data[0:8])[0],'gb2312') #站号
        pdata['siteName']=self.rhexgbk(struct.unpack("32s",block_data[8:40])[0],'gb2312') #站点名称
        pdata['latitude']=struct.unpack("f",block_data[40:44])[0] #纬度
        pdata['longitude']=struct.unpack("f",block_data[44:48])[0] #经度
        pdata['antennaHeight']=struct.unpack("i",block_data[48:52])[0] #天线高度
        pdata['groundHeight']=struct.unpack("i",block_data[52:56])[0] #地面高度
        pdata['frequency']=struct.unpack("f",block_data[56:60])[0] #工作频率
        pdata['beamWidthHori']=struct.unpack("f",block_data[60:64])[0] #水平波束宽度
        pdata['beamWidthVert']=struct.unpack("f",block_data[64:68])[0] #垂直波束宽度
        pdata['rdaVersion']=struct.unpack("i",block_data[68:72])[0] #RDA版本号
        pdata['radarType']=struct.unpack("h",block_data[72:74])[0] #雷达类型
        pdata['antennaGain']=struct.unpack("h",block_data[74:76])[0] #天线增益
        pdata['transmittingFeederLoss']=struct.unpack("h",block_data[76:78])[0] #发射馈线损耗
        pdata['receivingFeederLoss']=struct.unpack("h",block_data[78:80])[0] #接收馈线损耗
        pdata['otherLoss']=struct.unpack("h",block_data[80:82])[0] #其他损耗
        pdata['reserved']=self.rhexgbk(struct.unpack("46s",block_data[82:128])[0]) #保留字段
        self.data['common_block']['site_conf']=pdata
        # 任务配置信息
        block_data=f.read(256)
        pdata={}
        pdata['taskName']=self.rhexgbk(struct.unpack("32s",block_data[0:32])[0],'gb2312') #任务名称
        pdata['taskDescription']=self.rhexgbk(struct.unpack("128s",block_data[32:160])[0],'gb2312') #任务描述
        pdata['polarizationType']=struct.unpack("i",block_data[160:164])[0] #极化方式
        pdata['scanType']=struct.unpack("i",block_data[164:168])[0] #扫描任务类型
        pdata['pulseWidth']=struct.unpack("i",block_data[168:172])[0] #脉冲宽度
        pdata['scanStartTime']=struct.unpack("i",block_data[172:176])[0] #扫描开始时间
        pdata['cutNumber']=struct.unpack("i",block_data[176:180])[0] #扫描层数
        pdata['horizontalNoise']=struct.unpack("f",block_data[180:184])[0] #水平通道噪声
        pdata['verticalNoise']=struct.unpack("f",block_data[184:188])[0] #垂直通道噪声
        pdata['horizontalCalibration']=struct.unpack("f",block_data[188:192])[0] #水平通道标定值
        pdata['verticalCalibration']=struct.unpack("f",block_data[192:196])[0] #垂直通道标定值
        pdata['horizontalNoiseTemperature']=struct.unpack("f",block_data[196:200])[0] #水平通道噪声温度
        pdata['verticalNoiseTemperature']=struct.unpack("f",block_data[200:204])[0] #垂直通道噪声温度
        pdata['zdrCalibration']=struct.unpack("f",block_data[204:208])[0] #ZDR标定偏差
        pdata['phidpCalibration']=struct.unpack("f",block_data[208:212])[0] #差分相移标定偏差
        pdata['ldrCalibration']=struct.unpack("f",block_data[212:216])[0] #系统LDR标定偏差
        pdata['reserved']=self.rhexgbk(struct.unpack("40s",block_data[216:256])[0]) #保留字段
        self.data['common_block']['task_conf']=pdata
        # 扫描层配置信息
        for i in range(self.data['common_block']['task_conf']['cutNumber']):
            block_data=f.read(256)
            pdata={}
            pdata['processMode']=struct.unpack("i",block_data[0:4])[0] #处理模式
            pdata['waveForm']=struct.unpack("i",block_data[4:8])[0] #波形类别
            pdata['prf1']=struct.unpack("f",block_data[8:12])[0] #脉冲重复频率1
            pdata['prf2']=struct.unpack("f",block_data[12:16])[0] #脉冲重复频率2
            pdata['dealiasingMode']=struct.unpack("i",block_data[16:20])[0] #速度退模糊方法
            pdata['azimuth']=struct.unpack("f",block_data[20:24])[0] #方位角
            pdata['elevation']=struct.unpack("f",block_data[24:28])[0] #俯仰角
            pdata['startAngle']=struct.unpack("f",block_data[28:32])[0] #起始角度
            pdata['endAngle']=struct.unpack("f",block_data[32:36])[0] #结束角度
            pdata['angularResolution']=struct.unpack("f",block_data[36:40])[0] #角度分辨率
            pdata['scanSpeed']=struct.unpack("f",block_data[40:44])[0] #扫描速度
            pdata['logResolution']=struct.unpack("i",block_data[44:48])[0] #强度分辨率
            pdata['dopplerResolution']=struct.unpack("i",block_data[48:52])[0] #多普勒分辨率
            pdata['maximumRange']=struct.unpack("i",block_data[52:56])[0]##1 最大距离1
            pdata['maximumRange']=struct.unpack("i",block_data[56:60])[0]##2 最大距离2
            pdata['startRange']=struct.unpack("i",block_data[60:64])[0] #起始距离
            pdata['sample']=struct.unpack("i",block_data[64:68])[0]##1 采样个数1
            pdata['sample']=struct.unpack("i",block_data[68:72])[0]##2 采样个数2
            pdata['phaseMode']=struct.unpack("i",block_data[72:76])[0] #相位编码模式
            pdata['atmosphericLoss']=struct.unpack("f",block_data[76:80])[0] #大气衰减
            pdata['nyquistSpeed']=struct.unpack("f",block_data[80:84])[0] #最大不模糊速度
            pdata['momentsMask']=struct.unpack("l",block_data[84:92])[0] #数据类型掩码
            pdata['momentsSizeMask']=struct.unpack("l",block_data[92:100])[0] #数据大小掩码
            pdata['miscFilterMask']=struct.unpack("i",block_data[100:104])[0] #滤波设置掩码
            pdata['sqiThreshold']=struct.unpack("f",block_data[104:108])[0] #SQI门限
            pdata['sigThreshold']=struct.unpack("f",block_data[108:112])[0] #SIG门限
            pdata['csrThreshold']=struct.unpack("f",block_data[112:116])[0] #CSR门限
            pdata['logThreshold']=struct.unpack("f",block_data[116:120])[0] #LOG门限
            pdata['cpaThreshold']=struct.unpack("f",block_data[120:124])[0] #CPA门限
            pdata['pmiThreshold']=struct.unpack("f",block_data[124:128])[0] #PMI门限
            pdata['dplogThreshold']=struct.unpack("f",block_data[128:132])[0] #PMI门限
            pdata['thresholdsr']=struct.unpack("4s",block_data[132:136])[0] #阈值门限保留
            pdata['dBTMask']=struct.unpack("i",block_data[136:140])[0] #dBT质控掩码
            pdata['dBZMask']=struct.unpack("i",block_data[140:144])[0] #dBZ质控掩码
            pdata['velocityMask']=struct.unpack("i",block_data[144:148])[0] #速度质控掩码
            pdata['spectrumWidthMask']=struct.unpack("i",block_data[148:152])[0] #谱宽质控掩码
            pdata['dpMask']=struct.unpack("i",block_data[152:156])[0] #偏振量质控掩码
            pdata['maskReserved']=self.rhexgbk(struct.unpack("12s",block_data[156:168])[0]) #质控掩码保留位
            pdata['reserved1']=self.rhexgbk(struct.unpack("4s",block_data[168:172])[0]) #扫描同步标志
            pdata['direction']=struct.unpack("i",block_data[172:176])[0] #天线运行方向
            pdata['groundClutterClassifierType']=struct.unpack("h",block_data[176:178])[0] #地物杂波图类型
            pdata['groundClutterFilterType']=struct.unpack("h",block_data[178:180])[0] #地物滤波类型
            pdata['groundClutterFilterNotchWidth']=struct.unpack("h",block_data[180:182])[0] #地物滤波宽度
            pdata['groundClutterFilterWindow']=struct.unpack("h",block_data[182:184])[0] #滤波窗口类型
            pdata['reserved2']=self.rhexgbk(struct.unpack("72s",block_data[184:256])[0]) #保留字段
            self.data['common_block']['cut_conf'].append(pdata)
        # 径向数据
        is_last_radial=False
        while not is_last_radial:
            # 径向头
            block_data=f.read(64)
            pdata={}
            pdata['radialState']=struct.unpack("i",block_data[0:4])[0] #径向数据状态
            pdata['spotBlank']=struct.unpack("i",block_data[4:8])[0] #消隐标志
            pdata['sequenceNumber']=struct.unpack("i",block_data[8:12])[0] #序号
            pdata['radialNumber']=struct.unpack("i",block_data[12:16])[0] #径向数
            pdata['elevationNumber']=struct.unpack("i",block_data[16:20])[0] #仰角编号
            pdata['azimuth']=struct.unpack("f",block_data[20:24])[0] #方位角
            pdata['elevation']=struct.unpack("f",block_data[24:28])[0] #仰角
            pdata['seconds']=struct.unpack("i",block_data[28:32])[0] #秒
            pdata['microseconds']=struct.unpack("i",block_data[32:36])[0] #微秒
            pdata['lengthofdata']=struct.unpack("i",block_data[36:40])[0] #数据长度
            pdata['momentNumber']=struct.unpack("i",block_data[40:44])[0] #数据类别数量
            pdata['reserved1']=struct.unpack("h",block_data[44:46])[0] #保留字段
            pdata['horizontalEstimatedNoise']=struct.unpack("h",block_data[46:48])[0] #径向的水平估计噪声
            pdata['verticalEstimatedNoise']=struct.unpack("h",block_data[48:50])[0] #径向的垂直估计噪声
            pdata['zipType']=struct.unpack("c",block_data[50:51])[0] #压缩类型
            pdata['reserved2']=self.rhexgbk(struct.unpack("13s",block_data[51:64])[0]) #保留字段
            pdata['datahead']={}
            for ri in range(pdata['momentNumber']):
                # 径向数据头
                block_data=f.read(32)
                dh={}
                dh['dataType']=struct.unpack("i",block_data[0:4])[0] #数据类型
                dh['scale']=struct.unpack("i",block_data[4:8])[0] #比例
                dh['offset']=struct.unpack("i",block_data[8:12])[0] #偏移
                dh['binLength']=struct.unpack("h",block_data[12:14])[0] #库字节长度
                dh['flags']=struct.unpack("h",block_data[14:16])[0] #标志
                dh['length']=struct.unpack("i",block_data[16:20])[0] #长度
                dh['reserved']=self.rhexgbk(struct.unpack("12s",block_data[20:32])[0]) #保留字段
                # 径向数据计算
                radial_numpy=f.read(dh['length'])
                if dh['binLength']==1:
                    radial_numpy=np.frombuffer(radial_numpy, dtype='u1').astype(float)
                elif dh['binLength']==2:
                    radial_numpy=np.frombuffer(radial_numpy, dtype='u2').astype(float)
                else:
                    raise RadarError("雷达格式手册未说明有binLength不属于[1,2]的情况，请核实!")
                if dh['scale']>0:
                    radial_numpy=(radial_numpy-dh['offset'])/dh['scale']
                else:
                    radial_numpy=(radial_numpy-dh['offset'])*dh['scale']
                # 将数据存在data['radial']
                datakey=self.get_datatype_byid(dh['dataType'])
                if datakey not in self.data['radial'].keys():
                    self.data['radial'][datakey]={
                        'data':[],
                        'dataLength':[],
                        'radialState':[],
                        'spotBlank':[],
                        'sequenceNumber':[],
                        'radialNumber':[],
                        'elevationNumber':[],
                        'azimuth':[],
                        'elevation':[],
                        'seconds':[],
                        'microseconds':[],
                        'horizontalEstimatedNoise':[],
                        'verticalEstimatedNoise':[],
                        'ppiElevation':[],
                        'rhiAzimuth':[]
                    }
                self.data['radial'][datakey]['data'].append(radial_numpy)
                self.data['radial'][datakey]['dataLength'].append(len(radial_numpy))
                self.data['radial'][datakey]['radialState'].append(pdata['radialState'])
                self.data['radial'][datakey]['spotBlank'].append(pdata['spotBlank'])
                self.data['radial'][datakey]['sequenceNumber'].append(pdata['sequenceNumber'])
                self.data['radial'][datakey]['radialNumber'].append(pdata['radialNumber'])
                self.data['radial'][datakey]['elevationNumber'].append(pdata['elevationNumber'])
                self.data['radial'][datakey]['azimuth'].append(pdata['azimuth'])
                self.data['radial'][datakey]['elevation'].append(pdata['elevation'])
                self.data['radial'][datakey]['seconds'].append(pdata['seconds'])
                self.data['radial'][datakey]['microseconds'].append(pdata['microseconds'])
                self.data['radial'][datakey]['horizontalEstimatedNoise'].append(pdata['horizontalEstimatedNoise'])
                self.data['radial'][datakey]['verticalEstimatedNoise'].append(pdata['verticalEstimatedNoise'])
#                 if self.data['common_block']['task_conf']['scanType'] in [0,1]:
#                     ele=self.getStandardElevation(pdata['elevation'],self.data['common_block']['task_conf']['taskName'])
#                     self.data['radial'][datakey]['ppiElevation'].append(ele)
#                 else:
#                     self.data['radial'][datakey]['rhiAzimuth'].append(self.data['common_block']['cut_conf'][pdata['elevationNumber']-1]['azimuth'])
                # 径向数据头存储在此
                dh['dataid']=len(self.data['radial'][datakey]['data'])-1 #数据存储位置
                pdata['datahead'][datakey]=dh
                
            # 是否是最后一个径向数据，如果要退出循环了
            if pdata['radialState']==4:
                is_last_radial=True
                   
        return self.data
    
    def get_standard_elevation(self,elevation,vcp='VCP21D'):
        """
        获取标准仰角
        
        对于同一仰角PPI扫描，不同射线给出的仰角值不同，这里做标准化处理。
        例外：例如VCP21D第一层二层角度一样，因此不论给的是第一层的角度还是第二层角度，返回的layerid都是第一层，因此layerid仅供参考
        
        
        Parameters:
        elevation (float): 射线真实角度
        vcp (str): VCP模式
        
        Return:
        (layerid,elevation)
        """
        std_list=self.vcp_standard_elevation[vcp]
        position=np.argmin(np.abs(np.array(std_list)-elevation))
        return (position+1,std_list[position])
    
    def get_standard_elevation_byid(self,layerid,vcp='VCP21D'):
        """
        通过层id获取标准仰角
        
        对于同一仰角PPI扫描，不同射线给出的仰角值不同，这里做标准化处理。
        例外：例如VCP21D第一层二层角度一样，因此不论给的是第一层的角度还是第二层角度，返回的layerid都是第一层，因此layerid仅供参考
        
        
        Parameters:
        layerid (int): 层id
        vcp (str): VCP模式
        
        Return:
        elevation
        """
        return self.vcp_standard_elevation[vcp][layerid-1]
        
        