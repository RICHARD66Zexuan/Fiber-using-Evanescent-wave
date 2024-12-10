#simple model 
#minangle=75


# ----------------------------------------------------------------------------
# Set parameters
# ----------------------------------------------------------------------------

dpi=3000 #If you want the jpg form of figure
increament=0.0035 # The interval of angle

# ----------------------------------------------------------------------------
# import module
# ----------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.tri import TriAnalyzer, Triangulation, UniformTriRefiner
import math
import matplotlib
from matplotlib.font_manager import FontProperties

# ----------------------------------------------------------------------------
# import custom module (in Subfunction branch)
# ----------------------------------------------------------------------------

import function as f
import function_simple as f1
import get_point 
import Add_boundary

# ----------------------------------------------------------------------------
# import font 
# ----------------------------------------------------------------------------

#An example, if you want to get the figure's font you want. 
#font_path='/content/drive/MyDrive/software_TiO2/font/times.ttf'
#matplotlib.font_manager.fontManager.addfont(font_path)
#font_properties = FontProperties(fname="/content/drive/MyDrive/software_TiO2/font/calibri.ttf", size=10)

# ----------------------------------------------------------------------------
# Set file address
# ----------------------------------------------------------------------------

#The simulated experimental data you want. (In xlsx)
experimental_file='/content/drive/MyDrive/software_TiO2/Excel/Exper/nT/nTchange_simple.xlsx'

#Just the path of file is OK.
#First
workbook_file1='/content/drive/MyDrive/software_TiO2/Excel/Patchiness_Za/pza_RSM_169boundaru1.xlsx'
#Second
workbook_file2='/content/drive/MyDrive/software_TiO2/Excel/Patchiness_Za/pza_RSM_169boundaru2.xlsx'
#First
workbook_file_delaunay1='/content/drive/MyDrive/software_TiO2/Excel/result_RSM_169boundary_ttt.xlsx'
#Second
workbook_file_delaunay2='/content/drive/MyDrive/software_TiO2/Excel/result_RSM_169boundary_tttt.xlsx'


# ----------------------------------------------------------------------------
# Get two new Excel
# ----------------------------------------------------------------------------

def create_excel_xls(path):
    data_df = pd.DataFrame()
    writer = pd.ExcelWriter(path)
    data_df.to_excel(writer)
    writer.close()

create_excel_xls(workbook_file1)
writer1=pd.ExcelWriter(workbook_file1)

create_excel_xls(workbook_file2)
writer2=pd.ExcelWriter(workbook_file2)

# ----------------------------------------------------------------------------
# Delaunay parameters(First)
# ----------------------------------------------------------------------------

subdiv_1 =0
init_mask_frac_1 = 0
min_circle_ratio_1 = 0.04


# ----------------------------------------------------------------------------
# Delaunay parameters(Second)
# ----------------------------------------------------------------------------

subdiv_2 =0
init_mask_frac_2 = 0
min_circle_ratio_2 = 0.04


# ----------------------------------------------------------------------------
# Set getting points' parameter(First)
# ----------------------------------------------------------------------------

p_x1 =0.0001
za_y1 =0.000000001

interval1_x1=0.16
interval2_y1=0.000000200

n_x1= 15
m_y1= 40

# ----------------------------------------------------------------------------
# Set getting points' parameter(Second)
# ----------------------------------------------------------------------------

p_x2 =0.15
za_y2 =0.000000001

interval1_x2=0.85
interval2_y2=0.000000200

n_x2= 15
m_y2= 15

# ----------------------------------------------------------------------------
# Generate points (First)
# ----------------------------------------------------------------------------

interval_pacthiness = interval1_x1/n_x1
interval_za = interval2_y1/m_y1

point1 = []
point2 = []

for i in range(n_x1) :
    
    for j in range(m_y1):
    
        point11 = get_point.get_RSM_point1(p_x1 + i * interval_pacthiness , p_x1 + (i+1) * interval_pacthiness  )
        point1 = point1 + point11
        
        point22 = get_point.get_RSM_point2(za_y1 + j * interval_za , za_y1 + (j+1) * interval_za  )
        point2 = point2 + point22

boundary=Add_boundary.boundary_point(p_x1,p_x1+interval1_x1,za_y1,za_y1+interval2_y1,interval_pacthiness,interval_za)
point1=point1+boundary[0]+[p_x1+interval1_x1]
point2=point2+boundary[1]+[za_y1+interval2_y1]


df=pd.DataFrame({'patchiness':point1 , 'za':point2})

df.to_excel(writer1)
writer1.close()

# ----------------------------------------------------------------------------
# Generate points (Second)
# ----------------------------------------------------------------------------

interval_pacthiness = interval1_x2/n_x2
interval_za = interval2_y2/m_y2

point1 = []
point2 = []

for i in range(n_x2) :
    
    for j in range(m_y2):
    
        point11 = get_point.get_RSM_point1(p_x2 + i * interval_pacthiness , p_x2 + (i+1) * interval_pacthiness  )
        point1 = point1 + point11
        
        point22 = get_point.get_RSM_point2(za_y2 + j * interval_za , za_y2 + (j+1) * interval_za  )
        point2 = point2 + point22

boundary=Add_boundary.boundary_point(p_x2,p_x2+interval1_x2,za_y2,za_y2+interval2_y2,interval_pacthiness,interval_za)
point1=point1+boundary[0]+[p_x2+interval1_x2]
point2=point2+boundary[1]+[za_y2+interval2_y2]


df=pd.DataFrame({'patchiness':point1 , 'za':point2})

df.to_excel(writer2)
writer2.close()



# ----------------------------------------------------------------------------
# Get Result Excel for Delaunay (prepare parameters)
# ----------------------------------------------------------------------------

#Define exper variable
medium_refraction_index=[]
length=[]
wavelength=[]
diameter=[]
incident_energy=[]
fiber_index=[]
nT=[]
lamp_radius=[]
max_distance=[]
k=[]

df=pd.read_excel(experimental_file)
number=len(df.values) #The number of Exper
for i in range(number):
    medium_refraction_index.append(df.values[i,0])
    length.append(df.values[i,1])
    diameter.append(df.values[i,2])
    wavelength.append(df.values[i,3])
    incident_energy.append(df.values[i,4])
    fiber_index.append(df.values[i,5])
    nT.append(df.values[i,6])
    lamp_radius.append(df.values[i,7])
    max_distance.append(df.values[i,8])
    k.append(df.values[i,9])

df1=pd.read_excel(workbook_file1) 
df2=pd.read_excel(workbook_file2) 




# ----------------------------------------------------------------------------
# Get Result Excel for Delaunay (First)
# ----------------------------------------------------------------------------

p=[]
za=[]
create_excel_xls(workbook_file_delaunay1)
writer=pd.ExcelWriter(workbook_file_delaunay1)
m=len(df1.values)
for i in range(m): #读取数据
    p.append(df1.values[i,1])
    za.append(df1.values[i,2])
for i in [0]:  #主函数
    minangle=1.30899693
    maxangle=1.569
    fiber_critical_angle=f.get_fiber_critical_angle(fiber_index[i], medium_refraction_index[i])
    n=f.calculate_angle_divided_number(minangle, maxangle, increament)
    incident_energy[i]=f1.chageincident_energy_simple(incident_energy[i], n)
    angle=f.get_corresponding_angle(n, minangle, maxangle, fiber_critical_angle, increament)
    edeep=f.calculate_edeep(fiber_index[i], angle, medium_refraction_index[i], wavelength[i])
    Tqt=f.calculate_Tqt(n, fiber_index[i], angle, nT[i])
    Tmedium=f.calculate_Tmedium(medium_refraction_index[i], fiber_index[i], angle)
    N=f1.calculate_reflection_number_simple(n,angle, length[i], diameter[i])
    z=f.calculate_intermediate_quantity_z(m, angle, p, za, length[i], diameter[i], Tqt, incident_energy[i], edeep, N)
    energy_edisspated=f.calculate_energy_edisspated_list(m, angle, incident_energy[i], z)
    energy_rdisspated=f.calculate_energy_rdisspated_list(m, angle, incident_energy[i], z)
    energy_nontir=f.calculate_energy_nontir(m, angle, p, Tqt, incident_energy[i], N, Tmedium)
    disspated_energy=f1.sum_disspated_energy_123_simple(energy_edisspated, energy_rdisspated, m,energy_nontir)
    disspated_radiao=f1.calculate_radio_13_mode123_simple(m,n, incident_energy[i] , disspated_energy)
    radio_e=f1.calculate_radio_e_simple(energy_edisspated,m, incident_energy[i],disspated_energy)
    radio_r=f1.calculate_radio_r_simple(energy_rdisspated, m, incident_energy[i],disspated_energy)
    absorb_coefficient=f.calculate_absorb_coefficient(wavelength[i],k[i],0.0000001)
    near_field_ratio=[radio_e[k]+radio_r[k]*absorb_coefficient for k in range(m)]
    print('1down')
    
    df=pd.DataFrame({'patchiness':p,'za':za,'E_edis_ratio':radio_e[i]})
    df.to_excel(writer)
    writer.close()



# ----------------------------------------------------------------------------
# Get Result Excel for Delaunay (Second)
# ----------------------------------------------------------------------------

p=[]
za=[]
create_excel_xls(workbook_file_delaunay2)
writer=pd.ExcelWriter(workbook_file_delaunay2)
m=len(df2.values)
for i in range(m): #读取数据
    p.append(df2.values[i,1])
    za.append(df2.values[i,2])
for i in [0]:  #主函数
    minangle=1.30899693
    maxangle=1.569
    fiber_critical_angle=f.get_fiber_critical_angle(fiber_index[i], medium_refraction_index[i])
    n=f.calculate_angle_divided_number(minangle, maxangle, increament)
    incident_energy[i]=f1.chageincident_energy_simple(incident_energy[i], n)
    angle=f.get_corresponding_angle(n, minangle, maxangle, fiber_critical_angle, increament)
    edeep=f.calculate_edeep(fiber_index[i], angle, medium_refraction_index[i], wavelength[i])
    Tqt=f.calculate_Tqt(n, fiber_index[i], angle, nT[i])
    Tmedium=f.calculate_Tmedium(medium_refraction_index[i], fiber_index[i], angle)
    N=f1.calculate_reflection_number_simple(n,angle, length[i], diameter[i])
    z=f.calculate_intermediate_quantity_z(m, angle, p, za, length[i], diameter[i], Tqt, incident_energy[i], edeep, N)
    energy_edisspated=f.calculate_energy_edisspated_list(m, angle, incident_energy[i], z)
    energy_rdisspated=f.calculate_energy_rdisspated_list(m, angle, incident_energy[i], z)
    energy_nontir=f.calculate_energy_nontir(m, angle, p, Tqt, incident_energy[i], N, Tmedium)
    disspated_energy=f1.sum_disspated_energy_123_simple(energy_edisspated, energy_rdisspated, m,energy_nontir)
    disspated_radiao=f1.calculate_radio_13_mode123_simple(m,n, incident_energy[i] , disspated_energy)
    radio_e=f1.calculate_radio_e_simple(energy_edisspated,m, incident_energy[i],disspated_energy)
    radio_r=f1.calculate_radio_r_simple(energy_rdisspated, m, incident_energy[i],disspated_energy)
    absorb_coefficient=f.calculate_absorb_coefficient(wavelength[i],k[i],0.0000001)
    near_field_ratio=[radio_e[k]+radio_r[k]*absorb_coefficient for k in range(m)]
    print('2down')
    
    df=pd.DataFrame({'patchiness':p,'za':za,'E_edis_ratio':radio_e[i]})
    df.to_excel(writer)
    writer.close()



# ----------------------------------------------------------------------------
# Get Result Excel for Delaunay (First)
# ----------------------------------------------------------------------------

df=pd.read_excel(workbook_file_delaunay1)

# Initial points
x_initial=df['patchiness'].values
y_initial=df['za'].values
z_initial=df['E_edis_ratio'].values

x=np.array(x_initial)
y=np.array(y_initial)
z=np.array(z_initial)

# meshing with Delaunay triangulation
tri = Triangulation(x, y)
ntri = tri.triangles.shape[0]

# masking badly shaped triangles at the border of the triangular mesh.
mask = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio_1)
tri.set_mask(mask)

# refining the data
refiner = UniformTriRefiner(tri)
tri_refi, z_test_refi = refiner.refine_field(z, subdiv=subdiv_1)

# 输出新的xy点阵
p1=tri_refi.x
za1=tri_refi.y

m1=len(p1)  #patchiness和za的对数
print(m1)



# ----------------------------------------------------------------------------
# Get Result Excel for Delaunay (Second)
# ----------------------------------------------------------------------------

df=pd.read_excel(workbook_file_delaunay2)

# Initial points
x_initial=df['patchiness'].values
y_initial=df['za'].values
z_initial=df['E_edis_ratio'].values

x=np.array(x_initial)
y=np.array(y_initial)
z=np.array(z_initial)

# meshing with Delaunay triangulation
tri = Triangulation(x, y)
ntri = tri.triangles.shape[0]

# masking badly shaped triangles at the border of the triangular mesh.
mask = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio_2)
tri.set_mask(mask)

# refining the data
refiner = UniformTriRefiner(tri)
tri_refi, z_test_refi = refiner.refine_field(z, subdiv=subdiv_2)

# 输出新的xy点阵
p2=tri_refi.x
za2=tri_refi.y

m2=len(p2)  #patchiness和za的对数
print(m2)

# ----------------------------------------------------------------------------
# Get data
# ----------------------------------------------------------------------------

za11=za1*1000000000
za22=za2*1000000000
tri1=Triangulation(p1, za11)
tri2=Triangulation(p2, za22)
print(za11)
print(za1)



# ----------------------------------------------------------------------------
# Calculation(First)(Second)
# ----------------------------------------------------------------------------

for i in range(number):  #主函数
    minangle=1.30899693
    maxangle=1.569
    fiber_critical_angle=f.get_fiber_critical_angle(fiber_index[i], medium_refraction_index[i])
    n=f.calculate_angle_divided_number(minangle, maxangle, increament)
    incident_energy[i]=f1.chageincident_energy_simple(incident_energy[i], n)
    angle=f.get_corresponding_angle(n, minangle, maxangle, fiber_critical_angle, increament)
    edeep=f.calculate_edeep(fiber_index[i], angle, medium_refraction_index[i], wavelength[i])
    Tqt=f.calculate_Tqt(n, fiber_index[i], angle, nT[i])
    Tmedium=f.calculate_Tmedium(medium_refraction_index[i], fiber_index[i], angle)
    N=f1.calculate_reflection_number_simple(n,angle, length[i], diameter[i])
    absorb_coefficient=f.calculate_absorb_coefficient(wavelength[i],k[i],0.0000001)
    
    z=f.calculate_intermediate_quantity_z(m1, angle, p1, za1, length[i], diameter[i], Tqt, incident_energy[i], edeep, N)
    energy_edisspated=f.calculate_energy_edisspated_list(m1, angle, incident_energy[i], z)
    energy_rdisspated=f.calculate_energy_rdisspated_list(m1, angle, incident_energy[i], z)
    energy_nontir=f.calculate_energy_nontir(m1, angle, p1, Tqt, incident_energy[i], N, Tmedium)
    del z
    disspated_energy=f1.sum_disspated_energy_123_simple(energy_edisspated, energy_rdisspated, m1,energy_nontir)
    disspated_radiao1=f1.calculate_radio_13_mode123_simple(m1,n, incident_energy[i] , disspated_energy)
    radio_e1=f1.calculate_radio_e_simple(energy_edisspated,m1, incident_energy[i],disspated_energy)
    radio_r1=f1.calculate_radio_r_simple(energy_rdisspated, m1, incident_energy[i],disspated_energy)
    del energy_edisspated
    del energy_rdisspated
    del energy_nontir
    del disspated_energy
    
    near_field_ratio1=[radio_e1[k]+radio_r1[k]*absorb_coefficient for k in range(m1)]
    print('1第{number}组实验'.format(number=i+1))
    
    z=f.calculate_intermediate_quantity_z(m2, angle, p2, za2, length[i], diameter[i], Tqt, incident_energy[i], edeep, N)
    energy_edisspated=f.calculate_energy_edisspated_list(m2, angle, incident_energy[i], z)
    energy_rdisspated=f.calculate_energy_rdisspated_list(m2, angle, incident_energy[i], z)
    energy_nontir=f.calculate_energy_nontir(m2, angle, p2, Tqt, incident_energy[i], N, Tmedium)
    del z
    del N
    disspated_energy=f1.sum_disspated_energy_123_simple(energy_edisspated, energy_rdisspated, m2,energy_nontir)
    disspated_radiao2=f1.calculate_radio_13_mode123_simple(m2,n, incident_energy[i] , disspated_energy)
    radio_e2=f1.calculate_radio_e_simple(energy_edisspated,m2, incident_energy[i],disspated_energy)
    radio_r2=f1.calculate_radio_r_simple(energy_rdisspated, m2, incident_energy[i],disspated_energy)
    del energy_edisspated
    del energy_rdisspated
    del energy_nontir
    del disspated_energy
    
    near_field_ratio2=[radio_e2[k]+radio_r2[k]*absorb_coefficient for k in range(m2)]
    print('2第{number}组实验'.format(number=i+1))
    
    
    levels = np.arange(0, 1.1, 0.1)
    my_x_ticks=np.arange(0, 1.1, 0.1)
    # ----------------------------------------------------------------------------
    # Graph Z1
    # ----------------------------------------------------------------------------
    
    z1=disspated_radiao1
    del disspated_radiao1
    
    fig,ax=plt.subplots()
    
    cs=ax.tricontourf(tri1, z1, levels=levels,cmap='jet')
    #operate X
    plt.xlim(0,1)
    
    plt.xticks(my_x_ticks,fontsize=10, fontproperties=font_properties)
    plt.xlabel(r"$\mathit{p}$ (cm$^2$/cm$^2$)", fontsize=10, fontproperties=font_properties)
    plt.yticks(fontsize=10, fontproperties=font_properties)
    plt.ylabel(r"$\mathit{Z}$$_{a}$ (nm)",fontsize=10, fontproperties=font_properties)
    
    
    z1=disspated_radiao2
    del disspated_radiao2
    
    cs=ax.tricontourf(tri2, z1, levels=levels,cmap='jet')
    
    plt.title(f"$n_c$ = {nT[i]:.2f}", fontsize=10, fontproperties=font_properties)
    #set colorbar
    cbar=fig.colorbar(cs)
    colorbarticks=np.arange(0, 1.1, 0.1)
    cbar.set_ticks(colorbarticks,labels=['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'],fontsize=10, fontproperties=font_properties)
    cbar.ax.set_title(r"$\mathit{E}_{\mathrm{dis}}$/$\mathit{E}_{\mathrm{i}}$", fontsize=10, fontproperties=font_properties)
   
    
    
    # ----------------------------------------------------------------------------
    # Save z1
    # ----------------------------------------------------------------------------
    
    #plt.savefig('/content/drive/MyDrive/software_TiO2/Excel/Simple_Result/z1.svg',format='svg')
    #plt.savefig('/content/drive/MyDrive/software_TiO2/Excel/Simple_Result/z1.png')
    plt.savefig(f'/content/drive/MyDrive/software_TiO2/Result/nT_change/Z1/{i+1}.svg',format='svg')
    
    # ----------------------------------------------------------------------------
    # Graph Ze
    # ----------------------------------------------------------------------------
    
    ze=radio_e1
    del radio_e1
    
    fig,ax=plt.subplots()
    
    cs=ax.tricontourf(tri1, ze, levels=levels,cmap='jet')
    #operate X
    plt.xlim(0,1)
    
    plt.xticks(my_x_ticks,fontsize=10, fontproperties=font_properties)
    plt.xlabel(r"$\mathit{p}$ (cm$^2$/cm$^2$)", fontsize=10, fontproperties=font_properties)
    plt.yticks(fontsize=10, fontproperties=font_properties)
    plt.ylabel(r"$\mathit{Z}$$_{a}$ (nm)",fontsize=10, fontproperties=font_properties)
    
    ze=radio_e2
    del radio_e2
    
    cs=ax.tricontourf(tri2, ze, levels=levels,cmap='jet')
    
    plt.title(f"$n_c$ = {nT[i]:.2f}", fontsize=10, fontproperties=font_properties)
    
    #set colorbar
    cbar=fig.colorbar(cs)
    colorbarticks=np.arange(0, 1.1, 0.1)
    cbar.set_ticks(colorbarticks,labels=['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'],fontsize=10, fontproperties=font_properties)
    cbar.ax.set_title(r"$\mathit{E}_{\mathrm{Edis}}$/$\mathit{E}_{\mathrm{i}}$", fontsize=10, fontproperties=font_properties)
    
    
    
    # ----------------------------------------------------------------------------
    # Save ze
    # ----------------------------------------------------------------------------
    
    #plt.savefig('/content/drive/MyDrive/software_TiO2/Excel/Simple_Result/ze.png')
    #plt.savefig('/content/drive/MyDrive/software_TiO2/Excel/Simple_Result/ze.svg',format='svg')
    plt.savefig(f'/content/drive/MyDrive/software_TiO2/Result/nT_change/Ze/{i+1}.svg',format='svg')
    
    # ----------------------------------------------------------------------------
    # Graph ZNF
    # ----------------------------------------------------------------------------
    
    zNF=near_field_ratio1
    del near_field_ratio1
    
    fig,ax=plt.subplots()
    
    cs=ax.tricontourf(tri1, zNF, levels=levels,cmap='jet')
    #operate X
    plt.xlim(0,1)
    
    plt.xticks(my_x_ticks,fontsize=10, fontproperties=font_properties)
    plt.xlabel(r"$\mathit{p}$ (cm$^2$/cm$^2$)", fontsize=10, fontproperties=font_properties)
    plt.yticks(fontsize=10, fontproperties=font_properties)
    plt.ylabel(r"$\mathit{Z}$$_{a}$ (nm)",fontsize=10, fontproperties=font_properties)
    
    
    zNF=near_field_ratio2
    del near_field_ratio2
    
    cs=ax.tricontourf(tri2, zNF, levels=levels,cmap='jet')
    
    plt.title(f"$n_c$ = {nT[i]:.2f}", fontsize=10, fontproperties=font_properties)
    #set colorbar
    cbar=fig.colorbar(cs)
    colorbarticks=np.arange(0, 1.1, 0.1)
    cbar.set_ticks(colorbarticks,labels=['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'],fontsize=10, fontproperties=font_properties)
    cbar.ax.set_title(r"$\mathit{E}_{\mathrm{NF}}$/$\mathit{E}_{\mathrm{i}}$", fontsize=10, fontproperties=font_properties)
    
    
    
    # ----------------------------------------------------------------------------
    # Save zNF
    # ----------------------------------------------------------------------------

    plt.savefig(f'/content/drive/MyDrive/software_TiO2/Result/nT_change/ZNF/{i+1}.svg',format='svg')
    #plt.savefig('/content/drive/MyDrive/software_TiO2/Excel/Simple_Result/znf.png')
    #plt.savefig('/content/drive/MyDrive/software_TiO2/Excel/Simple_Result/znf.svg',format='svg')
