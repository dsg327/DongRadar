# -*- coding:utf-8 -*-
import numpy as np

def polar2cartesian(rho_list, theta_list, elevation_list, pole=[0, 0]):
    '''
    将极坐标系转换为卡迪尔坐标系

    https://wenku.baidu.com/view/bb2dc7d250e2524de5187e3e.html

    (极坐标单根射线坐标轴-极轴,极坐标弧度坐标轴,雷达层仰角,极点的经纬度坐标[latitude维度,longitude经度])    
    '''
    latx = list()
    lonx = list()
    height = list()
    earth_rm = 8500  # 等效地球半径
    for i,theta in enumerate(theta_list):
        for rho in rho_list:
            #print("angle:%s,distance:%s,elev:%s"%(theta,rho,elevation))
            x = np.cos(np.deg2rad(theta)) * rho * np.cos(np.deg2rad(elevation_list[i]))
            y = np.sin(np.deg2rad(theta)) * rho * np.cos(np.deg2rad(elevation_list[i]))
            deltalat = x / (111 * 1000)  # 球面弧度距离转角度值，111为1度对应的弧度
            lats = deltalat + pole[0]
            # deltalon = y / (111 * 1000)  # 直接用111则与天水雷达生成的图不一致（漳县回波偏东）
            deltalon = y / (111 * 1000 * np.cos(np.pi*lats/180)) #用维度修正距离之后绘制的图不再是原型，而是越往两极越扩大的类椭圆
            lons = deltalon + pole[1]
            h = rho * np.sin(np.deg2rad(elevation_list[i])) + rho ** 2 / (2 * earth_rm * 1000)  
            latx.append(lats)
            lonx.append(lons)
            height.append(h)
    return np.array(lonx), np.array(latx), np.array(height)


import math

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    使用Haversine公式计算两个经纬度点之间的球面距离。
    
    参数:
    lat1, lon1 : float
        第一个点的纬度和经度（十进制度数）。
    lat2, lon2 : float
        第二个点的纬度和经度（十进制度数）。
    
    返回:
    float
        两点之间的球面距离（公里）。
    """
    # 将经纬度转换为弧度
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # 计算纬度和经度的差值
    delta_lat = lat2 - lat1
    delta_lon = lon2 - lon1

    # 应用Haversine公式
    a = math.sin(delta_lat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(delta_lon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # 地球平均半径，单位公里
    R = 6371.0
    distance = R * c

    return distance

def calculate_bearing(lat1, lon1, lat2, lon2):
    """
    计算两个经纬度点之间的方位角。
    
    参数:
    lat1, lon1 : float
        第一个点的纬度和经度（十进制度数）。
    lat2, lon2 : float
        第二个点的纬度和经度（十进制度数）。
    
    返回:
    float
        两点之间的方位角（度）。
    """
    # 将经纬度转换为弧度
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # 计算经度和纬度的差值
    delta_lon = lon2 - lon1

    # 使用 atan2 来计算方位角
    bearing = math.degrees(math.atan2(math.sin(delta_lon), math.cos(lat1) * math.tan(lat2) - math.sin(lat1) * math.cos(delta_lon)))

    # 将方位角转换为 0-360 度的范围
    bearing = (bearing + 360) % 360

    return bearing