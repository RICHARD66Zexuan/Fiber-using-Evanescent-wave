import pandas as pd

import math

import numpy as np

from scipy.optimize import fsolve

def get_initial_critical_angle(fiber_index,medium_refraction_index,max_distance,lamp_radius):
    theta=[0.00000314159]
    c = lamp_radius - max_distance
    for i in range(999999):
        theta.append(theta[i]+0.00000314159*0.5)
    for i in range(1000000):
        a2 = math.asin(c * math.sin(theta[i]) / lamp_radius)
        if fiber_index * math.sin(a2) - medium_refraction_index > 0:
            break
    fiber_critical_angle=theta[i-1]
    return fiber_critical_angle

def get_maxa1(lamp_radius, max_distance, diameter, length, medium_refraction_index, fiber_index,fiber_critical_angle):
    divide = 100000  # 分的份数
    piece = (3.14159 * 0.5)/divide
    a1 = 3.14159/divide
    c = lamp_radius - max_distance
    halfpi = math.pi * 0.5
    b = 0
    for i in range(divide):
        a2 = math.asin(c * math.sin(a1) / lamp_radius)
        a3 = math.asin(math.sin(a2) * fiber_index / medium_refraction_index)
        a4 = a1 - a2
        a6 = halfpi - ((halfpi - (a1 - a2)) - a3)
        l1 = math.tan(a4) * lamp_radius
        l2 = (math.sin(a3) / math.cos(a6)) * ((l1 / math.sin(a4)) - lamp_radius)
        andl = l1 + l2
        a1 = a1 + piece
        if andl - diameter * 0.5 > 0:  # 这里这个标度如何还有待考虑
            b = a1
            break
    maxa1 = min(b - piece,fiber_critical_angle) #这里后续我们还需要去检验maxa1对应的minangle有没有超过了极小值点
    return maxa1

def angle_min(maxa1,lamp_radius, max_distance, diameter, length, medium_refraction_index, fiber_index):
    halfpi = math.pi * 0.5
    a2 = math.asin((lamp_radius - max_distance) * math.sin(maxa1) / lamp_radius) 
    a3 = math.asin(math.sin(a2) * fiber_index / medium_refraction_index)
    a4 = maxa1 - a2
    a6 = halfpi - ((halfpi - (a4)) - a3)
    a7 = math.asin(medium_refraction_index * math.sin(a6) / fiber_index)
    minangle = halfpi - a7
    #print("minangle=", end='')
    #print(minangle)
    return minangle
    
def get_mina1(lamp_radius, max_distance, diameter, length, medium_refraction_index, fiber_index,maxa1):
    divide = 100000  # 分的份数
    piece = (maxa1) / divide
    a1 = maxa1
    c = lamp_radius - max_distance
    halfpi = math.pi * 0.5
    for i in range(divide):
        a2 = math.asin(c * math.sin(a1) / lamp_radius)
        a3 = math.asin(math.sin(a2) * fiber_index / medium_refraction_index)
        a4 = a1 - a2
        a6 = halfpi - ((halfpi - (a1 - a2)) - a3)
        a7 = math.asin(medium_refraction_index * math.sin(a6) / fiber_index)
        l1 = math.tan(a4) * lamp_radius
        l2 = (math.sin(a3) / math.cos(a6)) * ((l1 / math.sin(a4)) - lamp_radius) 
        angle = halfpi - a7
        if math.tan(angle) >= (length / abs(0.5 * diameter - l1 - l2)):  #只要证明一点，maxa1可以保证不过图像的间断点,abs可以用
            break
        a1 = a1 - piece
    mina1=a1+piece
    return mina1

def get_maxangle(lamp_radius, max_distance, diameter, length, medium_refraction_index, fiber_index,mina1):
    halfpi = math.pi * 0.5
    c = lamp_radius - max_distance
    a2 = math.asin(c * math.sin(mina1) / lamp_radius)
    a3 = math.asin(math.sin(a2) * fiber_index / medium_refraction_index)
    a6 = halfpi - ((halfpi - (mina1 - a2)) - a3)
    a7 = math.asin(medium_refraction_index * math.sin(a6) / fiber_index)
    maxangle=halfpi-a7
    #print("maxangle=", end='')
    #print(maxangle)
    return maxangle

def get_fiber_critical_angle(fiber_index,medium_refraction_index):
    fiber_critical_angle=math.asin(medium_refraction_index/fiber_index)
    return fiber_critical_angle

def calculate_angle_divided_number(minangle, maxangle,increament):
    n = abs(int((maxangle - minangle)/increament))  # 预将最大最小角之间分成n份
    #print("total divided number=", end='')
    #print(n)
    return n  

def calculate_angle_divided_total_number(maxangle, increament):
    total_number_un= abs(int((0.5*math.pi-maxangle) / increament))
    #print("third group angles number=", end='')
    #print(total_number_un)
    return total_number_un

def chageincident_energy(incident_energy, n, total_number_un):
    incident_energy = incident_energy / (n + total_number_un)
    return incident_energy

def get_corresponding_angle(n, minangle, maxangle, fiber_critical_angle,increament):
    angle_initial=[minangle] #所有可能的角度
    angle1=[] #用于存储可以发生mode1与mode3的角度
    angle2=[] #用于存储可以发生mode1与mode2的角度
    for i in range(n - 1): #n-1保证了最大的angle小于等于maxangle
        angle_initial.append(angle_initial[i] + increament)
    for i in range(n):
        if angle_initial[i] < fiber_critical_angle:
            angle1.append(angle_initial[i])
        if angle_initial[i] >= fiber_critical_angle:
            angle2.append(angle_initial[i])
    angle=[]
    angle.append(angle1)
    angle.append(angle2)
    return angle #返回的是一个二维列表
    

def calculate_edeep(fiber_index, angle, medium_refraction_index, wavelength): #我们只对angle2中可以发生全反射的角求渗透深度
    edeep=[]
    n=len(angle[1])
    for i in range(n):
        m1 = math.sqrt(fiber_index * fiber_index * math.sin(angle[1][i]) * math.sin(angle[1][i]) - medium_refraction_index ** 2)
        edeep.append(wavelength / (4 * math.pi * m1))
    return edeep  

def calculate_Tqt(n, fiber_index, angle, nT):
    Tqt = []
    x=len(angle[0])
    for i in range(x):  
        a1=fiber_index * math.cos(angle[0][i])
        a2=nT * math.cos(angle[0][i])
        a3=1 - (fiber_index * fiber_index) * (math.sin(angle[0][i]) *math.sin(angle[0][i])) / (nT * nT)
        a4=math.sqrt(a3)
        b1=(a1 - nT * a4) / (a1 + nT * a4)
        b2=(fiber_index * a4 - a2) /(fiber_index * a4 +a2)
        Tqt.append(1 - 0.5 * (b1 * b1 + b2 * b2))    
    for i in range(n-x):  
        a1=fiber_index * math.cos(angle[1][i])
        a2=nT * math.cos(angle[1][i])
        a3=1 - (fiber_index * fiber_index) * (math.sin(angle[1][i]) *math.sin(angle[1][i])) / (nT * nT)
        a4=math.sqrt(a3)
        b1=(a1 - nT * a4) / (a1 + nT * a4)
        b2=(fiber_index * a4 - a2) /(fiber_index * a4 +a2)
        Tqt.append(1 - 0.5 * (b1 * b1 + b2 * b2))                     
    return Tqt

def calculate_Tmedium(medium_refraction_index,fiber_index,angle):
    Tmedium=[]
    n=len(angle[0])
    for i in range(n):  
        a1=fiber_index * math.cos(angle[0][i])
        a2=medium_refraction_index * math.cos(angle[0][i])
        a3=1 - (fiber_index * fiber_index) * (math.sin(angle[0][i]) *math.sin(angle[0][i])) / (medium_refraction_index * medium_refraction_index)
        a4=math.sqrt(a3)
        b1=(a1 - medium_refraction_index * a4) / (a1 + medium_refraction_index * a4)
        b2=(fiber_index * a4 - a2) /(fiber_index * a4 +a2)
        Tmedium.append(1 - 0.5 * (b1 * b1 + b2 * b2))
    return Tmedium

def f(x, theta, medium_refraction_index, fiber_index, max_distance ,lamp_radius):
    return np.array([abs(x[0]) + theta-math.pi / 2,
                     medium_refraction_index * math.sin(abs(x[1])) - fiber_index * math.sin(abs(x[0])),
                     abs(x[2]) + abs(x[3]) - math.pi/2,
                     abs(x[4]) - abs(x[1]) + abs(x[3]),
                     ((lamp_radius) / math.sin(abs(x[2]))) - (abs(x[5]) / math.sin(abs(x[3]))),
                     ((lamp_radius) / math.sin(abs(x[2]))) - lamp_radius - (abs(x[6])*math.sin(abs(x[1])) / math.sin(abs(x[4]))),
                     (fiber_index / medium_refraction_index) * (lamp_radius-max_distance)/math.sin(abs(x[4])) -
                     ((lamp_radius) /math.sin(abs(x[3]) + math.asin((medium_refraction_index / fiber_index)*math.sin(abs(x[4])))))])

def l1l2(n, angle, medium_refraction_index, fiber_index, lamp_radius, max_distance):
    list11=[]
    for i in range(n):
        angle_single=angle[i]
        result_fsolve= fsolve(f, [0.3, 0.4, 0.2, 0.2, 0.2, 0.0005, 0.0005], (angle_single, medium_refraction_index, fiber_index, max_distance, lamp_radius))
        list11.append(result_fsolve)
    return list11

def calculate_reflection_number(n,angle, length, diameter, fiber_index, medium_refraction_index, lamp_radius, max_distance):
    N=[]
    distance=[]
    x=len(angle[0]) #确定循环的间断点
    m1=l1l2(x, angle[0], medium_refraction_index, fiber_index, lamp_radius, max_distance)
    m2=l1l2(n-x, angle[1], medium_refraction_index, fiber_index, lamp_radius, max_distance)
    for i in range(x):
        distance.append(0.5 * diameter - abs(m1[i][5]) - abs(m1[i][6]))
    for i in range(n-x):
        distance.append(0.5 * diameter - abs(m2[i][5]) - abs(m2[i][6]))
    for i in range(x):
        N.append(int((((-1 * distance[i]) * math.tan(angle[0][i]) ) + length) / (diameter * math.tan(angle[0][i])))+1)
    for i in range(n-x):
        N.append(int((((-1 * distance[i]) * math.tan(angle[1][i]) ) + length) / (diameter * math.tan(angle[1][i])))+1)
    #print(N)
    return N

def calculate_intermediate_quantity_z(m, angle, p, za, length, diameter, Tqt, incident_energy, edeep, N):
    x=len(angle[0])
    n=len(angle[1])
    z=[]
    z1=[]
    z2=[]
    z3=[]
    z4=[]
    z11=[]
    z22=[]
    z33=[]
    z44=[]
    for i in range(m): #做p和za的循环
        for j in range(n): #做全反射相应角个数的循环
            a1= math.exp(-1 * (za[i] / edeep[j]))
            b1= (1 - p[i]) * a1
            z33.append(b1)
            b2= p[i] * Tqt[x+j]
            z44.append(b2)
            c1= ((1 - p[i]) * (1 - a1) + p[i] - b2) ** N[x+j] #应该是x+j
            z11.append(c1)
            c2= b1 + b2
            z22.append(c2)
        z1.append(z11)
        z2.append(z22)
        z3.append(z33)
        z4.append(z44)
        z11=[]
        z22=[]
        z33=[]
        z44=[]
    z.append(z1)
    z.append(z2)
    z.append(z3)
    z.append(z4)
    return z

def calculate_energy_edisspated_list(m, angle, incident_energy,z):  
    n=len(angle[1])
    energy_edisspated=[]
    for i in range(m): #做p和za的大循环
        temporary_sum=0
        for j in range(n): #做全反射个数的小循环
             temporary_sum=temporary_sum + (z[2][i][j] * (1 -z[0][i][j])) / z[1][i][j]
        energy_edisspated.append(incident_energy * temporary_sum)
    return energy_edisspated

def calculate_energy_rdisspated_list(m, angle, incident_energy,z):
    n=len(angle[1])
    energy_rdisspated=[]
    for i in range(m): #做p和za的大循环
        temporary_sum=0
        for j in range(n): #做全反射个数的小循环
             temporary_sum=temporary_sum + (z[3][i][j] * (1 -z[0][i][j])) / z[1][i][j]
        energy_rdisspated.append(incident_energy * temporary_sum) 
    return energy_rdisspated
             
def calculate_energy_nontir(m, angle, p, Tqt, incident_energy, N ,Tmedium):
    n=len(angle[0])
    #print("first group angle number=", end='')
    #print(n)
    energy_nontir=[]
    for i in range(m): #做p和za的大循环
        temporary_sum=0
        for j in range(n): #做非全反射个数的小循环
            temporary_sum=temporary_sum + (1 - ((1 - p[i]) * (1 - Tmedium[j]) + p[i] * (1 - Tqt[j])) ** N[j] )
        energy_nontir.append(incident_energy * temporary_sum)
    return energy_nontir


#不知道有什么用地函数奇奇怪怪地
def calculate_energy_initial_output(n, angle,incident_energy, N ,Tmedium, total_number_un):
    c=len(angle[0])
    energy_initial_output=incident_energy * (n + total_number_un)
    a1=[]
    d=0
    for i in range(c):
        f = 0
        a1.append(incident_energy * Tmedium[i])
        for j in range(N[i]):
           if j == 0 :
               a = a1[i]
           else: a = (incident_energy - f) * Tmedium[i]
           f = f + a
        d= d + f   
    energy_initial_output=energy_initial_output - d 
    return energy_initial_output

def sum_disspated_energy_13(energy_edisspated, energy_rdisspated, m):
    disspated_energy=[]
    for i in range(m): #做p和za的大循环
        sum1=energy_edisspated[i] + energy_rdisspated[i]
        disspated_energy.append(sum1)
    return disspated_energy

def sum_disspated_energy_123(energy_edisspated, energy_rdisspated, m,energy_nontir):
    disspated_energy=[]
    for i in range(m): #做p和za的大循环
        sum1=energy_edisspated[i] + energy_rdisspated[i] + energy_nontir[i]
        disspated_energy.append(sum1)
    return disspated_energy

def calculate_radio_2(m,angle,  energy_initial_output, energy_nontir, incident_energy):
    disspated_radio=[]
    n = len(angle[0])
    c = incident_energy * n
    for i in range(m):
        x= energy_nontir[i] / c
        disspated_radio.append(x)
    return  disspated_radio
        
def calculate_radio_13(m, total_number_un,angle, incident_energy , disspated_energy):
    disspated_radio=[]
    c= len(angle[1])
    #print("second and third group angles number=", end='')
    #print(c)
    for i in range(m):
        x=disspated_energy[i] / (incident_energy * (c + total_number_un) )
        disspated_radio.append(x)
    return disspated_radio

#专门用于处理Calculation_Mode1_Angle.py的情况
def calculate_radio_13_mode123(m, total_number_un,angle, incident_energy , disspated_energy):
    disspated_radio=[]
    a= len(angle[0])
    c= len(angle[1])
    #print("second and third group angles number=", end='')
    print(a)
    print(total_number_un)
    print(c)
    for i in range(m):
        x=disspated_energy[i] / (incident_energy * (a+c + total_number_un) )
        disspated_radio.append(x)
    return disspated_radio

def  calculate_radio_e(energy_edisspated,m, total_number_un,angle, incident_energy,disspated_energy):

    radio_e=[]
    c= len(angle[1])
    for i in range(m):
        x=energy_edisspated[i] / disspated_energy[i]
        radio_e.append(x)
    return radio_e

def  calculate_radio_r(energy_rdisspated, m, total_number_un,angle, incident_energy,disspated_energy):

    radio_r=[]
    c= len(angle[1])
    for i in range(m):
        x=energy_rdisspated[i] / disspated_energy[i]
        radio_r.append(x)
    return radio_r

def  calculate_absorb_coefficient(wavelength,k,s):

    coefficient=math.exp(-4 * math.pi * k * s / wavelength)
    c= 1- coefficient
    return c

