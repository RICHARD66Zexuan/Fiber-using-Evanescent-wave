import numpy as np

#Fill the boundary of the mesh created by Response surface methodology

def boundary_point(A_initial,A_final,B_initial,B_final,A_interval,B_interval):
  point1=np.arange(A_initial,A_final,A_interval).tolist()
  point2=np.arange(B_initial,B_final,B_interval).tolist()
  n1=len(point1)
  n2=len(point2)
  point1_add1=n2*[A_initial]
  point2_add1=n1*[B_initial]
  point1_add2=n2*[A_final]
  point2_add2=n1*[B_final]
  point=[]
  p=point1
  p=p+point1_add1
  p=p+point1
  p=p+point1_add2
  
  Za=point2_add1
  Za=Za+point2
  Za=Za+point2_add2
  Za=Za+point2

  point.append(p)
  point.append(Za)
  return point
