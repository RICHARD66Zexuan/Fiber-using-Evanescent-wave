import math
def chageincident_energy_simple(incident_energy, n):
    
    incident_energy = incident_energy / n
    return incident_energy

def calculate_reflection_number_simple(n,angle, length, diameter):
    
    N=[]
    x=len(angle[0]) 

    for i in range(x):
        N.append(int(length/(diameter * math.tan(angle[0][i])))+1)
    for i in range(n-x):
        N.append(int(length/(diameter * math.tan(angle[1][i])))+1)
    return N

def calculate_radio_13_mode123_simple(m,n, incident_energy , disspated_energy):
    
    disspated_radio=[]
    
    for i in range(m):
        x=disspated_energy[i] / (incident_energy * n )
        disspated_radio.append(x)
    return disspated_radio

def sum_disspated_energy_123_simple(energy_edisspated, energy_rdisspated, m,energy_nontir):
    
    disspated_energy=[]
    
    for i in range(m): 
        sum1=energy_edisspated[i] + energy_rdisspated[i] + energy_nontir[i]
        disspated_energy.append(sum1)
    return disspated_energy

def  calculate_radio_e_simple(energy_edisspated,m, incident_energy,disspated_energy):

    radio_e=[]
    
    for i in range(m):
        x=energy_edisspated[i] / disspated_energy[i]
        radio_e.append(x)
    return radio_e

def  calculate_radio_r_simple(energy_rdisspated, m, incident_energy,disspated_energy):

    radio_r=[]
    
    for i in range(m):
        x=energy_rdisspated[i] / disspated_energy[i]
        radio_r.append(x)
    return radio_r
