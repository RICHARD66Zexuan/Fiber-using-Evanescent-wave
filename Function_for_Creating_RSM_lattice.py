import math

def get_RSM_point1(high1, low1):
    average1 = 0.5 * (high1 + low1)

    a= 0.25 * (2 - math.sqrt(2))
    
    b= 0.25 * (2 + math.sqrt(2))
    
    a1 = a * low1 + b * high1
    
    b1 = b * low1 + a * high1
    
    point1 =[]
    
    point1.append(average1)
    
    point1.append(b1)
    
    point1.append(a1)
    
    point1.append(average1)
    
    point1.append(high1)
      
    point1.append(b1)
     
    point1.append(a1)
     
    point1.append(average1)
      
    point1.append(low1)
    
    return point1

def get_RSM_point2(high2, low2):
    average2 = 0.5 * (high2 + low2)
    
    a= 0.25 * (2 - math.sqrt(2))
    
    b= 0.25 * (2 + math.sqrt(2))

    a2 = a * low2 + b * high2

    b2 = b * low2 + a * high2
    
    point2 = []
    
    point2.append(average2)
    
    point2.append(b2)
    
    point2.append(b2)
    
    point2.append(high2)
    
    point2.append(average2)
    
    point2.append(a2)
    
    point2.append(a2)
    
    point2.append(low2)
    
    point2.append(average2)
    
    return point2

