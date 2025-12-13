import numpy as np
import mpmath as mp
import random as rand
import datetime
import time
import os
import subprocess
import matplotlib.pyplot as plt
from array import array

from ROOT import TMath, TF1, gRandom, TRandom

def dis_probability(theta, I_0):
    return I_0 * np.cos(theta)

### ------------- Funciones de distribución angular theta ----------- ##
def dis_thet(lx, lpar): ## Distribucion dis_angular
    theta = lx[0]
    return (TMath.Cos(theta)**2) * (TMath.Sin(theta))

def random_thet():
    gRandom.SetSeed(0); ## Cambia la semilla aleatoria para el GetRandom
    f = TF1("", dis_thet, 0, np.pi/2,1)
    #f = TF1("", "cos(x)^2 * sin(x)", 0, np.pi/2, 1)
    # print(type(f))
    thet = f.GetRandom()
    return thet 
### ----------------------------------------------------------------- ###

def coord_cartesian(Thet, Phi): ##### Coordenadas Cartesianas
    coord_X = np.sin(Thet) * np.cos(Phi)
    coord_Y = np.sin(Thet) * np.sin(Phi)
    coord_Z = np.cos(Thet)

    Vec_sph = [float(coord_X), float(coord_Y), float(coord_Z)]
    return Vec_sph 

def norma_vec(Vector):
    norma = np.sqrt(Vector[0] ** 2 + Vector[1] ** 2  + Vector[2] ** 2 )
    return norma

def Energy_list_log(lim_inf, lim_sup, N):
    list_Energy = []
    Delta_log = (np.log10(lim_sup) - np.log10(lim_inf))/ (N - 1)

    for i in np.arange(0, N):
        En = 10 ** (np.log10(lim_inf) + i * Delta_log)
        list_Energy.append(En)

    return list_Energy

def dimension_x(long_x, Step):
    step = Step

    list_long_x = [-long_x]

    while long_x:
        x = np.around(list_long_x[-1] + step, 4)
        list_long_x.append(x)
        # print(x)
        if list_long_x[-1] == long_x:
            break

    return list_long_x

def dimension_y(long_y, Step):
    step = Step

    list_long_y = [-long_y]

    while long_y:
        y = np.around(list_long_y[-1] + step, 4)
        # print(y)
        list_long_y.append(y)

        if list_long_y[-1] == long_y:
            break

    return list_long_y

def dimension_z(long_z, Step):
    step = Step

    list_long_z = [-long_z]
    # list_long_z = [0]
    # print(long_z)

    while long_z:
        z = np.around(list_long_z[-1] + step, 5)
        # print(z)
        list_long_z.append(z)

        if list_long_z[-1] == long_z:
            break

    return list_long_z

def intersection_CCD(flags_CCD, list_z, medida_z, Random_th ):
    # list_delta_L = []
    # list_z = [z_3, z_4, z_5, z_6]

    delta_L = None

    ## Caras 1 y 2
    if flags_CCD[0] and flags_CCD[1]:
        delta_L = (2 * medida_z) / np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        n_muons_in_CCD = 1
        # list_delta_L.append(delta_L)
        # continue

    ## Caras 1 y 3
    if flags_CCD[0] and flags_CCD[2]:
        h = medida_z - list_z[0]
        delta_L = h /  np.cos(Random_th)  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 4
    if flags_CCD[0] and flags_CCD[3]:
        h = medida_z - list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 5
    if flags_CCD[0] and flags_CCD[4]:
        h = medida_z - list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 1 y 6
    if flags_CCD[0] and flags_CCD[5]:
        h = medida_z - list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 2
    if flags_CCD[2] and flags_CCD[1]:
        h = medida_z + list_z[0] 
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 4
    if (list_z[0] > list_z[1]) and flags_CCD[2] and flags_CCD[3]:
        h = list_z[0] - list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 5
    if (list_z[0] > list_z[2]) and flags_CCD[2] and flags_CCD[4]:
        h = list_z[0] - list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 3 y 6
    if (list_z[0] > list_z[3]) and flags_CCD[2] and flags_CCD[5]:
        h = list_z[0] - list_z[3]
        delta_L =  (h/  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 2
    if flags_CCD[3] and flags_CCD[1]:
        h = medida_z + list_z[1] 
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 3
    if (list_z[1] > list_z[0]) and flags_CCD[3] and flags_CCD[2]:
        h = list_z[1] - list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 5
    if (list_z[1] > list_z[2]) and flags_CCD[3] and flags_CCD[4]:
        h = list_z[1] - list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 4 y 6
    if (list_z[1] > list_z[3]) and flags_CCD[3] and flags_CCD[5]:
        h = list_z[1] - list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 2
    if flags_CCD[4] and flags_CCD[1]:
        h = list_z[2] + medida_z
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 3
    if (list_z[2] > list_z[0]) and flags_CCD[4] and flags_CCD[2]:
        h = list_z[2] - list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 4
    if (list_z[2] > list_z[1]) and flags_CCD[4] and flags_CCD[3]:
        h = list_z[2] - list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 5 y 6
    if (list_z[2] > list_z[3]) and flags_CCD[4] and flags_CCD[5]:
        h = list_z[2]  - list_z[3]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 2
    if flags_CCD[5] and flags_CCD[1]:
        h = list_z[3] + medida_z
        delta_L =  (h /np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 3
    if (list_z[3] > list_z[0]) and flags_CCD[5] and flags_CCD[2]:
        h = list_z[3] - list_z[0]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 4
    if (list_z[3] > list_z[1]) and flags_CCD[5] and flags_CCD[3]:
        h = list_z[3] - list_z[1]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    ## Caras 6 y 5
    if (list_z[3] > list_z[2]) and flags_CCD[5] and flags_CCD[4]:
        h = list_z[3] - list_z[2]
        delta_L =  (h /  np.cos(Random_th))  ## cm
        # n_muons_in_CCD = n_muons_in_CCD + 1
        # list_delta_L.append(delta_L)
        n_muons_in_CCD = 1
        # continue

    # if (flags_CCD[0] == False) and  (flags_CCD[1] == False) and (flags_CCD[2] == False) and (flags_CCD[3] == False) and (flags_CCD[4] == False) and (flags_CCD[5] == False):
    if delta_L is None:    
        delta_L, n_muons_in_CCD = 0, 0

    # if delta_L < 0.0725 and delta_L != 0:
    #     print(delta_L)

    # return np.abs(delta_L), n_muons_in_CCD
    return delta_L, n_muons_in_CCD

### ======================= Funcinones Smith-Duller (ROOT) ========================== ###
def Smith_Dull(lx , lpar): ### Modelo de Smith-Duller
    E_mu = lx[0]
    theta = lpar[0]

    ## ------------- Constantes físicas -------------- ##
    k = 8 / 3
    b = 0.771
    lambda_pi = 120     ## g/cm^2
    y_0 = 1000      ## g/cm^2
    r = 0.76
    b_mu = 0.8      

    # # if units == 0:
    a = 2.5         ## MeV cm^2/g
    m_mu = 105.64    ## MeV/c^2 
    m_pi = 139.6    ## MeV/c^2

    # elif units == 1:
    # a = 0.0025         ## GeV cm^2/g
    # m_mu = 0.1057    ## GeV/c^2 
    # m_pi = 0.1396    ## GeV/c^2

    tau_mu_0 = 2.2 * 10**(-6)   ## s
    tau_0 = 2.6 * 10 **(-8)     ## s
    rho_0 = 0.00129 ## g/cm^3
    c = TMath.C() * 100 ## cm/s

    ### ---------------------- Parámetros ---------------------- ###
    E_pi = (1 / r) * (E_mu + a * y_0 * ((1/TMath.Cos(theta)) - 0.1))
    B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c)
    # B_mu = (b_mu * m_mu * y_0 * c)/(tau_mu_0 * rho_0)
    P_mu = ((0.1 * TMath.Cos(theta)) * (1 - (a * (y_0 *(1/TMath.Cos(theta)) - 100))/( r * E_pi)) ) ** ((B_mu)/((r * E_pi + 100 * a) * TMath.Cos(theta)))
    j_pi = (m_pi * y_0)/(c * tau_0 * rho_0)
    # j_pi = (m_pi * y_0 * c)/(tau_0 * rho_0)

    ## Intensidad diferencial
    # C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    C_2 = E_pi * TMath.Cos(theta)
    C_3 = b * j_pi

    # return (C_1 * np.sin(theta)) / (C_2 + C_3)
    En_k = (C_1 ) / (C_2 + C_3)
    return En_k

def random_SD(theta):
    gRandom.SetSeed(0) ## Cambia la semilla aleatoria para el GetRandom
    f = TF1("", Smith_Dull, 1, 10**5,1)
    # f = TF1("", LandV, 0, 10, 2)
    f.SetParameter(0, theta)
    Enpri = f.GetRandom()
    return Enpri 
### ================================================================================= ###

### ========================= Funciones de Landau (ROOT) ============================ ###
def LandV(lx, lpar):
    Delta = lx[0]	# Energy loss in absorber
    L = lpar[0]		# Thickness of absorber (Distance crossed by the particle)
    p = lpar[1]		# Momentum (in MeV/c)

    K = 0.1535	# K coefficient = 2*pi*N*r^2*m*c^2 (in MeV g^-1 cm^2)

    z = -1;		# Charge number of incident particle
    ZA = 0.498487	# Atomic number over Atomic mass of absorber (for Si)
    c = TMath.C()	# Speed of light
    me = 0.510998928	#  Electron mass in MeV/c^2
    M = 105.65839	# Muon mass in MeV/c^2
    I = 0.000173		# Mean excitation energy in MeV (for Si)
    
    bg = p/M
    beta = bg/np.sqrt(1+bg**2)	# Beta factor
    gamma = 1/np.sqrt(1-beta**2)	# Gamma factor
    pi = TMath.Pi()
    rho = 2.33	# Density of material in gr/cm^3 (for Si)

    a = 0.1492	# Parameters (taken from reference of PDF for SI)
    k = 3.2546		
    X0 = 0.2014		
    X1 = 2.8715		
    C = 4.4351		
    d0 = 0.14;	# For Si 
    # d0 = 0 # Like Bryan

    X = np.log10(bg)    
    if X>=X1:
        d = 2 * np.log(10.0) * X - C

    elif X0<=X and X<X1:
        d = 2 * np.log(10.0) * X - C + a * ((X1-X)**k)
        
    elif X<X0:
        d = d0 * 10**(2*(X-X0))

    WM = 2 * me * ((beta*gamma)**2)/(1+(2*me*gamma/M)+(me/M)**2)   #Maximum energy tranfer (Bryan from PDG)
    # WM = (2 * me * (beta*gamma)**2)/(1+(2*me/M)*np.sqrt(1 + (beta*gamma)**2) + (me/M)**2)   #Maximum energy tranfer

    loge = np.log((1-beta**2)*(I**2)/(2*me*(beta**2))) + (beta**2) # log epsilon variable (with electron mass like Bryan)
    # loge = np.log((1-beta**2)*(I**2)/(2*M*(beta**2))) + (beta**2) # log epsilon variable (with muon mass)

    EC = 0.577	# Euler's constant

    ##### ===================== With Leo version ===============0 ###
    xi = (K)*rho*ZA*L*(z/beta)**2		# Xi variable 
    DeltaAv = xi * (np.log((2*me*(bg**2)*WM)/(I**2))-(2*beta**2)-(d)) # Mean energy loss (Bethe-Bloch from PDG)
    # DeltaAv = xi * (np.log(2*me*(gamma**2)*(p**2)*WM/(I**2))-(2*beta**2)-(d)) # Mean energy loss (Bethe-Bloch from Leo)
    ### ========================================================= ###

    Lambda = (Delta-xi*(np.log(xi)-loge+1-EC))/xi # Lambda parameter

    # Deltamp = xi*(np.log(xi/np.exp(loge))+0.198-d)		# Most probable energy loss (Leo)
    # Deltamp = xi*(np.log(xi) - loge + 0.198 - d)		# Most probable energy loss (Leo)
    # lambdamp = (Deltamp-xi*(np.log(xi)-loge+1-EC))/xi 

    Deltamp = xi * (np.log((2 * me * (bg**2))/I) + np.log(xi/I) + 0.2 - beta**2 - d)		# Most probable energy loss from PDG
    lambdamp = (Deltamp-xi*(np.log(xi)-loge+1-EC))/xi

    kappa = xi/WM		# Kappa ratio
    beta2 = beta**2
    sigma2 = (xi**2)*(1-beta2/2)/kappa		# Standard deviation for relativistic particles

    if kappa<=0.01:
        phi = TMath.Landau(Lambda, lambdamp, 1.0)
        return phi/xi

    elif 0.01<kappa and kappa<10:
        # vav = TMath.Vavilov(Delta-Deltamp, kappa, beta2)
        vav = TMath.Vavilov(Lambda - lambdamp, kappa, beta2)
        return vav

    else:
        gauss = TMath.Gaus(Delta, DeltaAv, np.sqrt(sigma2))
        return gauss

def random_LV(s, p):
    gRandom.SetSeed(0) ## Cambia la semilla aleatoria para el GetRandom
    f = TF1("", LandV, 0, 4, 2)
    # f = TF1("", LandV, 0, 10, 2)
    # f.SetParameter(0, np.double(s))
    # f.SetParameter(1, np.double(p))
    f.SetParameters(s, p)

    Edep = f.GetRandom()


    ## =============== Cálculo de kappa ================= ##
    K = 0.1535	# K coefficient = 2*pi*N*r^2*m*c^2 (in MeV g^-1 cm^2)

    z = -1;		# Charge number of incident particle
    ZA = 0.498487	# Atomic number over Atomic mass of absorber (for Si)
    c = TMath.C()	# Speed of light
    me = 0.510998928	#  Electron mass in MeV/c^2
    M = 105.65839	# Muon mass in MeV/c^2
    I = 0.000173		# Mean excitation energy in MeV (for Si)
    bg = p/M
    beta = bg/np.sqrt(1+bg**2)	# Beta factor
    gamma = 1/np.sqrt(1-beta**2)	# Gamma factor
    pi = TMath.Pi()
    rho = 2.33	# Density of material in gr/cm^3 (for Si)

    xi = (K)*rho*ZA* s *(z/beta)**2		# Xi variable
    WM = 2 * me * ((bg)**2)/(1+(2*me*gamma/M)+(me/M)**2)   #Maximum energy tranfer (Bryan from PDG)
    kappa = xi/WM		# Kappa ratio

    return Edep * 1000, kappa    #### Energía en KeV
### ================================================================================= ###

def muon_generator_1(number_thet):
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []

    nmuons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    m_mu = 105.7    ## MeV/c^2

    for i in np.arange(0,number_thet):
        gRandom.SetSeed(0)
        try:
            ### ================== Seleccion aleatoria de theta, phi y en_pri(Smith-Duller) =============== ###
            Random_th = random_thet() ## Escoje un ángulo segun la distribución de Theta_true en radianes
            Random_phi = rand.random()*(2 * np.pi) ## Escoje un ángulo phi uniforme en radianes
            Random_energy = random_SD(Random_th) ## Escoje una energía cinética segun la distribución de Smith-Duller en MeV
            ### =========================================================================================== ###

            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)

            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '                  \r')

        except:
            print('Hubo un fallo')
            continue



    # n_muons_in_CCD = len(list_delta_L)
    dict_muons =  {'Theta(Rad)': list_thet_in_CCD, 'Phi(Rad)' : list_phi_in_CCD, 
                   'Energy-SD(MeV)' : list_energy_pri_in_CCD} 

    return dict_muons

def muon_generator_2(number_thet, Radio, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z, half_plane_size):
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []
    list_nmuons = []
    list_delta_L = []
    list_energy_Landau = []

    nmuons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    m_mu = 105.7    ## MeV/c^2

    for i in np.arange(0,number_thet):
        gRandom.SetSeed(0)
        ### ================== Seleccion aleatoria de theta, phi y en_pri(Smith-Duller) =============== ###
        Random_th = random_thet() ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.random()* (2 * np.pi) ## Escoje un ángulo phi uniforme en radianes
        # print(Random_phi)
        Random_energy = random_SD(Random_th) ## Escoje una energía cinética segun la distribución de Smith-Duller en MeV
        ### =========================================================================================== ###

        ### ==================== Momento del muon ======================= ###
        En_tot = Random_energy + m_mu  
        momentum = np.sqrt(En_tot**2 - m_mu**2)
        ### ============================================================= ###
        
        ### =========== Vector de direccion del muon y punto al azar sobre el plano =========== ###
        # flag_z = True
        # while flag_z:
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.

        # normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)
        normal_Vec =  [-1 * np.sin(Random_th) * np.cos(Random_phi), 
                    -1 * np.sin(Random_th) * np.sin(Random_phi), 
                    -1 * np.cos(Random_th)] ## Es un vector normal unitario apuntando  hacia el centro de coordenadas
        # print('Normal_vec: ', normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), 
                    np.cos(Random_th) * np.sin(Random_phi),
                    -np.sin(Random_th)]
        
        vec_phi = [-np.sin(Random_phi), 
                   np.cos(Random_phi), 
                   0]

        random_a = -half_plane_size + rand.random() * 2 * half_plane_size ## Selecciona un valor uniforme para el parámetro a
        random_b = -half_plane_size + rand.random() * 2 * half_plane_size ##      ''      ''      ''      ''          ''    b

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]

        random_plane_point = [Point[0] + P_vector[0], 
                              Point[1] + P_vector[1], 
                              Point[2] + P_vector[2]]

            # if random_plane_point[2] > 0:
            #     flag_z = False
        

        # print('Punto sobre el plano: ', random_plane_point)
        ### ===================================================================================== ###


        ### ==================== Intersecciones con cada cara =====================  ####
        flag_cara_1, flag_cara_2, flag_cara_3 = False, False, False
        flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False 

        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 
        point_sup = [x_1, y_1, medida_z]

        #### Cara Inferior ###
        t_2 = (-medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2
        point_inf = [x_2, y_2, 0]

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3
        pointx_1 = [medida_x, y_3, z_3]

        ### Cara 4 ###
        t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4
        pointx_2 = [-medida_x, y_4, z_4]

        #### Caras en Y ###
        ### Cara 5 ###
        t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5
        pointy_1 = [x_5, medida_y, z_5]

        ### Cara 6 ###
        t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_6
        x_6 = random_plane_point[0] + normal_Vec[0] * t_6
        pointy_2 = [x_6, -medida_y, z_6]

        list_z = [z_3, z_4, z_5, z_6]
        # print(list_z)
        limit_around = 4
        n_flags = 0
        
        flag_faces = True
        while flag_faces:
            if np.around(x_1, limit_around) in mapeo_x and np.around(y_1, limit_around) in mapeo_y:
                flag_cara_1 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_2, limit_around) in mapeo_x and np.round(y_2, limit_around) in mapeo_y:
                flag_cara_2 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3, limit_around) in mapeo_y and np.around(z_3, limit_around + 1) in mapeo_z:
                flag_cara_3 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(y_4, limit_around) in mapeo_y and np.round(z_4, limit_around + 1) in mapeo_z:
                flag_cara_4 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)

            if np.around(x_5, limit_around) in mapeo_x and np.around(z_5, limit_around + 1) in mapeo_z:
                flag_cara_5 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_6, limit_around) in mapeo_x and np.round(z_6, limit_around + 1) in mapeo_z:
                flag_cara_6 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)
            else:
                flag_faces = False


        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]

        # print(list_flags, Random_th)
        if n_flags == 1:
            print(list_flags)
            print("x1: ", np.around(x_1, limit_around))
            print("y1: ", np.around(y_1, limit_around))
            print("x2: ", np.around(x_2, limit_around))
            print("y2: ", np.around(y_2, limit_around))
            print("y3: ", np.around(y_3, limit_around))
            print("z3: ", np.around(z_3, limit_around + 1))
            print("y4: ", np.around(y_4, limit_around))
            print("z4: ", np.around(z_4, limit_around + 1))
            print("x5: ", np.around(x_5, limit_around))
            print("z5: ", np.around(z_5, limit_around + 1))
            print("x6: ", np.around(x_6, limit_around))
            print("z6: ", np.around(z_6, limit_around + 1))
        graf_flag = False
        if graf_flag:

            ### =============== Dibuja el plano tangente sobre la esfera ============== ###
            ax = plt.figure().add_subplot(projection='3d')

            # P_plane = [-half_plane_size * vec_thet[0] + -half_plane_size  * vec_phi[0], 
            #         -half_plane_size * vec_thet[1] + -half_plane_size * vec_phi[1], 
            #         -half_plane_size * vec_thet[2] + -half_plane_size * vec_phi[2]]
            # Point_plt = [Point[0] + P_plane[0], Point[1] - P_plane[1], Point[2] + P_plane[2]]
            # ax.scatter(Point_plt[0], Point_plt[1], Point_plt[2], color='k') 

            # P_plane = [half_plane_size * vec_thet[0] + half_plane_size  * vec_phi[0], 
            #         half_plane_size * vec_thet[1] + half_plane_size * vec_phi[1], 
            #         half_plane_size * vec_thet[2] + half_plane_size * vec_phi[2]]
            # Point_plt = [Point[0] + P_plane[0], Point[1] + P_plane[1], Point[2] - P_plane[2]]
            # ax.scatter(Point_plt[0], Point_plt[1], Point_plt[2], color='k') 

            # P_plane = [half_plane_size * vec_thet[0] + half_plane_size  * vec_phi[0], 
            #         half_plane_size * vec_thet[1] + half_plane_size * vec_phi[1], 
            #         half_plane_size * vec_thet[2] + half_plane_size * vec_phi[2]]
            # Point_plt = [Point[0] + P_plane[0], Point[1] + P_plane[1], Point[2] + P_plane[2]]
            # ax.scatter(Point_plt[0], Point_plt[1], Point_plt[2], color='k') 

            # P_plane = [-half_plane_size * vec_thet[0] + -half_plane_size  * vec_phi[0], 
            #         -half_plane_size * vec_thet[1] + -half_plane_size * vec_phi[1], 
            #         -half_plane_size * vec_thet[2] + -half_plane_size * vec_phi[2]]
            # Point_plt = [Point[0] + P_plane[0], Point[1] + P_plane[1], Point[2] - P_plane[2]]
            # ax.scatter(Point_plt[0], Point_plt[1], Point_plt[2], color='k') 

            ### ================== Grafico d ela CCD =================== ###
            ### =========== Cara inferior =============== ###
            X1 = [-medida_x, medida_x, medida_x, -medida_x, -medida_x, medida_x]
            Y1 = [-medida_y, -medida_y, medida_y, medida_y, -medida_y, medida_y]
            Z1 = [0, 0, 0, 0, 0, 0]
            ax.plot(xs=X1, ys=Y1, zs=Z1, color='k')

            ### =========== Cara superior =============== ###
            X1 = [-medida_x, medida_x, medida_x, -medida_x, -medida_x]
            Y1 = [-medida_y, -medida_y, medida_y, medida_y, -medida_y]
            Z1 = [medida_z, medida_z, medida_z, medida_z, medida_z]
            ax.plot(xs=X1, ys=Y1, zs=Z1, color='k', label = 'Thet: ' + str(Random_th) + '\n Phi: ' + str(Random_phi))

            ### ========== Uniones ==================== ###
            x1= [-medida_x, -medida_x]
            y1 = [-medida_y, -medida_y]
            z1 = [0, medida_z]
            ax.plot(xs=x1, ys=y1, zs=z1, color='k')

            x1= [medida_x, medida_x]
            y1 = [-medida_y, -medida_y]
            z1 = [0, medida_z]
            ax.plot(xs=x1, ys=y1, zs=z1, color='k')

            x1= [medida_x, medida_x]
            y1 = [medida_y, medida_y]
            z1 = [0, medida_z]
            ax.plot(xs=x1, ys=y1, zs=z1, color='k')

            x1= [-medida_x, -medida_x]
            y1 = [medida_y, medida_y]
            z1 = [0, medida_z]
            ax.plot(xs=x1, ys=y1, zs=z1, color='k')

            ### ============= Proyecion del domo ================ ###
            # theta = np.linspace(0, 4 * np.pi, 100)
            # # z = np.linspace(-2, 2, 100)
            # z = 0
            # r = Radio
            # x = r * np.sin(theta)
            # y = r * np.cos(theta)
            # ax.plot(x, y, z, label='parametric curve', linestyle='dashed')
            # ax.scatter([0,0], [0,0], [0,Radio])


            # ax.scatter([Point[0]], [Point[1]], [Point[2]], color='b') ## Punto sobre la esfera
            # ax.plot([0, Point[0]], [0, Point[1]], [0, Point[2]], label = 'Thet: ' + str(Random_th) + '\n Phi: ' + str(Random_phi),
            #         color='r', linestyle='dashed') ## Recta del origien al punto sobre la esf.       

            # ax.scatter(random_plane_point[0], random_plane_point[1], random_plane_point[2], 
            #             color='orange') ## Punto random del plano
            scale = 1
            X = [random_plane_point[0]/scale, point_sup[0]]
            Y = [random_plane_point[1]/scale, point_sup[1]]
            Z = [random_plane_point[2]/scale, point_sup[2]]
            ax.plot(X, Y, Z, color='k')
            ax.scatter(point_sup[0], point_sup[1], point_sup[2], color = 'k')

            X = [random_plane_point[0]/scale, point_inf[0]]
            Y = [random_plane_point[1]/scale, point_inf[1]]
            Z = [random_plane_point[2]/scale, point_inf[2]]
            ax.plot(X, Y, Z, color='r')
            ax.scatter(point_inf[0], point_inf[1], point_inf[2], color = 'k')

            #pointx_1
            X = [random_plane_point[0]/scale, pointx_1[0]]
            Y = [random_plane_point[1]/scale, pointx_1[1]]
            Z = [random_plane_point[2]/scale, pointx_1[2]]
            ax.plot(X, Y, Z, color='r')
            ax.scatter(pointx_1[0], pointx_1[1], pointx_1[2], color = 'k')

            X = [random_plane_point[0]/scale, pointx_2[0]]
            Y = [random_plane_point[1]/scale, pointx_2[1]]
            Z = [random_plane_point[2]/scale, pointx_2[2]]
            ax.plot(X, Y, Z, color='r')
            ax.scatter(pointx_2[0], pointx_2[1], pointx_2[2], color = 'k')

            X = [random_plane_point[0]/scale, pointy_1[0]]
            Y = [random_plane_point[1]/scale, pointy_1[1]]
            Z = [random_plane_point[2]/scale, pointy_1[2]]
            ax.plot(X, Y, Z, color='r')
            ax.scatter(pointy_1[0], pointy_1[1], pointy_1[2], color = 'k')

            X = [random_plane_point[0]/scale, pointy_2[0]]
            Y = [random_plane_point[1]/scale, pointy_2[1]]
            Z = [random_plane_point[2]/scale, pointy_2[2]]
            ax.plot(X, Y, Z, color='r')
            ax.scatter(pointy_2[0], pointy_2[1], pointy_2[2], color = 'k')

            # X = [random_plane_point[0], x_1]
            # Y = [random_plane_point[1], y_1]
            # Z = [random_plane_point[2], medida_z]
            # ax.plot(X, Y, Z, color='g')


            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.legend()
            plt.show()
            ### ======================================================================== ###

        Delta_L, _ = intersection_CCD(list_flags, list_z, medida_z, Random_th)
        # print(Delta_L)
        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para Delta L: ', Fin-In)

        if Delta_L > 0:
            # print(list_flags)
            # print(list_z)
            # print(Delta_L)
            nmuons_in_CCD =  nmuons_in_CCD + 1
            ### ================== Calculo de la energía de Landau ================ ###
            # Random_energy_Landau = random_LV(s=Delta_L, p = momentum)
            ### =================================================================== ###

            # list_nmuons.append(n_muon)
            # list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)
            
            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '\r')

            # else:
            #     n_negative_long = n_negative_long + 1
            #     continue

        else:
            # print(Delta_L)
            Random_energy_Landau = 0 # En KeV
            Delta_L = 0

            # list_nmuons.append(n_muon)
            # list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)

            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '\r')



    dict_muons =  {'NMuon': list_nmuons, 'Theta(Rad)': list_thet_in_CCD, 'Phi(Rad)' : list_phi_in_CCD, 
                   'Energy-SD(MeV)' : list_energy_pri_in_CCD, 'Delta_L(cm)' : list_delta_L, 
                   'Energy_Landau(KeV)' : list_energy_Landau} 

    return dict_muons, nmuons_in_CCD 

def muon_generator_3(number_thet, Radio, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z, half_plane_size):
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []
    list_nmuons = []
    list_delta_L = []
    list_energy_Landau = []
    list_kappa = []

    nmuons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    max_long = np.sqrt((medida_x * 2)**2 + (medida_y * 2)**2 + (medida_z * 2)**2)
    m_mu = 105.7    ## MeV/c^2

    for i in np.arange(0,number_thet):
        # found = False
        # while not found:
        #     rand_y = rand.random()
        #     rand_x = rand.random() * (np.pi/2)
        #     if rand_y < dis_angular(rand_x):
        #         Random_th = rand_x
        #         found = True

        gRandom.SetSeed(0)
        ### ================== Seleccion aleatoria de theta, phi y en_pri(Smith-Duller) =============== ###
        Random_th = random_thet() ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.random()*(2 * np.pi) ## Escoje un ángulo phi uniforme en radianes
        Random_energy = random_SD(Random_th) ## Escoje una energía cinética segun la distribución de Smith-Duller en MeV
        ### =========================================================================================== ###

        ### ==================== Momento del muon ======================= ###
        En_tot = Random_energy + m_mu  
        momentum = np.sqrt(En_tot**2 - m_mu**2)
        ### ============================================================= ###
        
        ### =========== Vector de direccion del muon y punto al azar sobre el plano =========== ###
        # flag_z = True
        # while flag_z:
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.

        # normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)

        normal_Vec =  [-1 * np.sin(Random_th) * np.cos(Random_phi), 
                       -1 * np.sin(Random_th) * np.sin(Random_phi), 
                       -1 * np.cos(Random_th)] ## Es un vector normal unitario apuntando  hacia el centro de coordenadas

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), 
                    np.cos(Random_th) * np.sin(Random_phi),
                    -np.sin(Random_th)]
        
        vec_phi = [-np.sin(Random_phi), 
                   np.cos(Random_phi), 
                   0]

        random_a = -half_plane_size + rand.random() * 2 * half_plane_size ## Selecciona un valor uniforme para el parámetro a
        random_b = -half_plane_size + rand.random() * 2 * half_plane_size ##      ''      ''      ''      ''          ''    b

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]

        random_plane_point = [Point[0] + P_vector[0], 
                              Point[1] + P_vector[1], 
                              Point[2] + P_vector[2]]

        # print('Punto sobre el plano: ', random_plane_point)
        ### ===================================================================================== ###

        ### ==================== Intersecciones con cada cara =====================  ####
        flag_cara_1, flag_cara_2, flag_cara_3 = False, False, False
        flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False 

        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

        #### Cara Inferior ###
        t_2 = (-medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3

        ### Cara 4 ###
        t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4

        #### Caras en Y ###
        ### Cara 3 ###
        t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5

        ### Cara 4 ###
        t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_6
        x_6 = random_plane_point[0] + normal_Vec[0] * t_6

        
        list_z = [z_3, z_4, z_5, z_6]
        # print(list_z)
        limit_around = 4
        n_flags = 0
        
        flag_faces = True
        while flag_faces:
            if np.around(x_1, limit_around) in mapeo_x and np.around(y_1, limit_around) in mapeo_y:
                flag_cara_1 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.around(x_2, limit_around) in mapeo_x and np.around(y_2, limit_around) in mapeo_y:
                flag_cara_2 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3, limit_around) in mapeo_y and np.around(z_3, limit_around + 1) in mapeo_z:
                flag_cara_3 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.around(y_4, limit_around) in mapeo_y and np.around(z_4, limit_around + 1) in mapeo_z:
                flag_cara_4 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)

            if np.around(x_5, limit_around) in mapeo_x and np.around(z_5, limit_around + 1) in mapeo_z:
                flag_cara_5 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.around(x_6, limit_around) in mapeo_x and np.around(z_6, limit_around + 1) in mapeo_z:
                flag_cara_6 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)
            
            else:
                flag_faces = False
        
        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]

        # print(list_flags, Random_th)
        if n_flags == 1:
            print(list_flags)
            print("x1: ", np.around(x_1, limit_around))
            print("y1: ", np.around(y_1, limit_around))
            print("x2: ", np.around(x_2, limit_around))
            print("y2: ", np.around(y_2, limit_around))
            print("y3: ", np.around(y_3, limit_around))
            print("z3: ", np.around(z_3, limit_around + 1))
            print("y4: ", np.around(y_4, limit_around))
            print("z4: ", np.around(z_4, limit_around + 1))
            print("x5: ", np.around(x_5, limit_around))
            print("z5: ", np.around(z_5, limit_around + 1))
            print("x6: ", np.around(x_6, limit_around))
            print("z6: ", np.around(z_6, limit_around + 1))

        
        try:
            Delta_L, _ = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        except:
            print('Hubo un error con el cálculo de L')
            print(Delta_L)
            # continue
            break

        #print(Delta_L)
        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para Delta L: ', Fin-In)

        if 0 < Delta_L and Delta_L < max_long:
            # print(list_flags)
            nmuons_in_CCD =  nmuons_in_CCD + 1
            # print('Delta L: ', Delta_L)
            # print('Energy S-D: ', Random_energy)
            # print('Momentum: ', momentum)

            ### ================== Calculo de la energía de Landau ================ ###
            try:
                Random_energy_Landau, kappa = random_LV(s=Delta_L, p = momentum)
                # print('Energy Landau: ', Random_energy_Landau)

            except: 
                print('Hubo en error con el ćalculo de Landau')
                print(Random_energy_Landau, kappa)
                break
                # contnue
            ### =================================================================== ###

            # list_nmuons.append(n_muon)
            list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)
            list_kappa.append(kappa)
            
            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '               \r')

            # else:
            #     n_negative_long = n_negative_long + 1
            #     continue

        else:
            Random_energy_Landau = 0 # En KeV
            Delta_L = 0
            kappa = -10

            # list_nmuons.append(n_muon)
            list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)
            list_kappa.append(kappa)

            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '                  \r')



    # n_muons_in_CCD = len(list_delta_L)
    dict_muons =  {'NMuon': list_nmuons, 'Theta(Rad)': list_thet_in_CCD, 'Phi(Rad)' : list_phi_in_CCD, 
                   'Energy-SD(MeV)' : list_energy_pri_in_CCD, 'Delta_L(cm)' : list_delta_L, 
                   'Energy_Landau(KeV)' : list_energy_Landau, 'kappa' : list_kappa} 

    return dict_muons, nmuons_in_CCD 

def muon_generator_CLUSTER(Energy, number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, 
                  long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z):
    
    list_rand_thet = []
    list_rand_phi = []
    list_rand_thet_deg = []
    list_rand_phi_deg = []
    list_random_energy = []

    ## Another lists ###
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []

    list_P_vector = []
    list_random_point = []
    list_delta_L = []
    list_energy_Landau = []

    m_mu = 105.7

    n_muons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    for i in np.arange(0,number_thet):
        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi en radianes
        # print(Random_th[0])
        Random_th_deg = np.degrees(Random_th[0]) ## El ángulo theta en grados
        Random_phi_deg = np.degrees(Random_phi) ## El ángulo phi en grados

        # list_rand_thet.append(Random_th[0])
        # list_rand_phi.append(Random_phi)
        # list_rand_thet_deg.append(Random_th_deg)
        # list_rand_phi_deg.append(Random_phi_deg)

        # list_dis_Energy = []
        # for energy in Energy:   ## Aquí se crea la distribución de Smith-Duller en MeV
        #     dis_Energy = dis_energy(energy, Random_th[0])
        #     list_dis_Energy.append(dis_Energy)
        
        # list_dis_Energy = dis_energy(Energy, Random_th[0], units=1)

        # Random_energy = rand.choices(Energy, list_dis_Energy)[0] * 1000 ## Escoje una energía segun la distribución de Smith-Duller en MeV
        # list_random_energy.append(Random_energy[0])

        list_rand_thet_deg.append(Random_th_deg)
        list_rand_phi_deg.append(Random_phi_deg)
        # list_random_energy.append(Random_energy)

        ### Momento del muon ###
        # momentum = np.sqrt(Random_energy[0]**2 - m_mu**2)
        # momentum = Random_energy

        # os.environ["EN_SMITH"] = str(momentum)
        
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
                                                                                            ##  hacia el centro de coordenadas
        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
        # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
        # list_random_phi.append(Random_phi)


        random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
        random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

        # list_random_a.append(random_a)
        # list_random_b.append(random_b)

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]
        
        # list_P_vector.append(P_vector)

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
        # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

        # print(random_plane_point)
        # list_random_point.append(random_plane_point)

        #### Intersecciones con cada cara   ####
        
        #### Cara Superior ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

        #### Cara Inferior ###
        t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3

        ### Cara 4 ###
        t_4 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4

        #### Caras en Y ###
        ### Cara 3 ###
        t_5 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5

        ### Cara 4 ###
        t_6 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_4 
        x_6 = random_plane_point[0] + normal_Vec[0] * t_4

        list_z = [z_3, z_4, z_5, z_6]
        
        if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
            flag_cara_1 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
            flag_cara_2 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
            flag_cara_3 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
            flag_cara_4 = True
            # print('Bandera 2: ', flag_cara_2)

        if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
            flag_cara_5 = True
            # print('Bandera 1: ', flag_cara_1) 

        if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
            flag_cara_6 = True
            # print('Bandera 2: ', flag_cara_2)

        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
        
        Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        if Delta_L != 0:

            if Delta_L > 0 and Delta_L < 2.1:
                # list_delta_L.append(Delta_L)
                # print('Estoy agregando el deltaL')

                n_muons_in_CCD = n_muons_in_CCD + muon


                # os.environ["DELTA_L"] = str(Delta_L)

                # print(os.environ["DELTA_L"])
                # print('Entrando al proceso de ROOT')
                
                # new_env = subprocess.run(["root", "-l", "-b", "/home/icn/mausanram/Software/CodigosICN/Simulacion_ab_initio/LandauVavilov_Mau.C", "-q"], 
                #                              capture_output=True)

                # print('Energía de SMith-Duller: ', os.getenv("EN_SMITH"))
                # print(new_env.stdout)
                # print(new_env.stderr)

                # print(os.getenv('PATH'))
                # subprocess.run()
                # print(new_env.stdout)
                # Random_energy_Landau = float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1]) # En KeV
                # print(Random_energy_Landau)

                # print(float(new_env.stdout.decode('ascii').split('=')[-1].split(' ')[1]))

                # Random_energy_Landau = random_LV(s = Delta_L, p = momentum) # En KeV


                list_rand_thet.append(Random_th[0])
                list_rand_phi.append(Random_phi)
                list_rand_thet_deg.append(Random_th_deg)
                list_rand_phi_deg.append(Random_phi_deg)
                # list_random_energy.append(Random_energy)
                list_delta_L.append(Delta_L)
                # list_energy_Landau.append(Random_energy_Landau)
                
                muon_in_bucle += 1
                # print("El valor de EDEP", str(os.getenv('USERNAME')))

                # print(muon_in_bucle, Random_th, Random_phi, Random_energy[0], Random_energy_Landau)

                print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet * number_points_per_angle), end = '\r')

            else:
                n_negative_long = n_negative_long + 1
                continue

        else:
                continue
        
        # print(os.environ)

        
    Len_Landau = len(list_energy_Landau)
    Len_Theta = len(list_rand_thet)

    # if Len

    # dict_muons =  {'Random_Thet': list_rand_thet, 'Random_Phi' : list_rand_phi, 'Random_Energy' : list_random_energy, 'DeltaL' : list_delta_L} 

    dict_muons =  {'Theta(Rad)': list_rand_thet, 'Theta(Deg)': list_rand_thet_deg, 
                   'Phi(Rad)' : list_rand_phi, 'Phi(Deg)' : list_rand_phi_deg, 
                   'Energy-SD(MeV)' : list_random_energy, 'Delta_L' : list_delta_L,
                   'Energy_Landau' : list_energy_Landau} 

    return dict_muons, n_muons_in_CCD, n_negative_long

def func_longitud(number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, 
                  long_a, long_b, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z):
    
    list_P_vector = []
    list_random_point = []
    list_delta_L = []

    list_random_th = []
    list_random_phi = []

    n_muons_in_CCD = 0
    n_negative_long = 0

    for i in np.arange(0,number_thet):
        list_random_energy_per_plane = []
        list_points_per_plane = []

        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        # print('Random Thet: ', np.degrees(Random_th[0]))
        # print('Random Phi: ', np.degrees(Random_phi))
        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))

        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        # normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma) 
        normal_Vec =  [-1 * np.sin(Random_th) * np.cos(Random_phi), -1 * np.sin(Random_th) * np.sin(Random_phi), -1 * np.cos(Random_th)]  ## Es un vector normal unitario apuntando hacia el centro de coordenadas
        # print('Direction Vec: ', normal_Vec)

        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), -np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        for i in np.arange(0,number_points_per_angle):
            flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
            # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
            # list_random_phi.append(Random_phi)


            random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            # list_random_a.append(random_a)
            # list_random_b.append(random_b)

            P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                        random_a * vec_thet[1] + random_b * vec_phi[1], 
                        random_a * vec_thet[2] + random_b * vec_phi[2]]
            
            list_P_vector.append(P_vector)

            random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

            # print(random_plane_point)
            list_random_point.append(random_plane_point)

            #### Intersecciones con cada cara   ####
            # print('Coord. del punto sobre el plano: (', random_plane_point[0], random_plane_point[1], random_plane_point[2],')')

            #### Cara Superior ###

            t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
            x_1 = random_plane_point[0] + normal_Vec[0] * t_1[0] 
            y_1 = random_plane_point[1] + normal_Vec[1] * t_1[0] 

            #### Cara Inferior ###
            t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
            x_2 = random_plane_point[0] + normal_Vec[0] * t_2[0] 
            y_2 = random_plane_point[1] + normal_Vec[1] * t_2[0]

            ### Caras en X ###
            ### Cara 3 ###
            t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0]
            z_3 = random_plane_point[2] + normal_Vec[2] * t_3[0] 
            y_3 = random_plane_point[1] + normal_Vec[1] * t_3[0]

            ### Cara 4 ###
            t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
            z_4 = random_plane_point[2] + normal_Vec[2] * t_4[0] 
            y_4 = random_plane_point[1] + normal_Vec[1] * t_4[0]

            #### Caras en Y ###
            ### Cara 5 ###
            t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1]
            z_5 = random_plane_point[2] + normal_Vec[2] * t_5[0] 
            x_5 = random_plane_point[0] + normal_Vec[0] * t_5[0]

            ### Cara 4 ###
            t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
            z_6 = random_plane_point[2] + normal_Vec[2] * t_6[0] 
            x_6 = random_plane_point[0] + normal_Vec[0] * t_6[0]

            # print(random_plane_point[0][0])
            # x, y, z = np.indices((10, 10, 10))

            # draw cuboids in the top left and bottom right corners, and a link between
            # them
            # cube1 = (x > 0) & (y > 0) & (z > 0)
            # cube2 = (x >= 5) & (y >= 5) & (z >= 5)
            # link = abs(x - y) + abs(y - z) + abs(z - x) <= 2

            # combine the objects into a single boolean array
            # voxelarray = cube1 | cube2 | link
            # voxelarray = cube1 | cube2 
            # voxelarray = cube1

            # set the colors of each object
            # colors = np.empty(voxelarray.shape, dtype=object)
            # colors[link] = 'red'
            # colors[cube1] = 'blue'
            # colors[cube2] = 'green'

            # and plot everything
            # ax = plt.figure().add_subplot(projection='3d')
            # # ax.voxels(filled = voxelarray, facecolors=colors)
            # VecStart_x = [-50, -50, -50, -50]
            # VecStart_y = [-5, -5, 5, 5]
            # VecStart_z = [0, 10, 0, 10]
            # VecEnd_x = [50, 50, 50, 50]
            # VecEnd_y = [-5, -5, 5, 5]
            # VecEnd_z  =[0, 10, 0, 10]

            # for i in range(0, len(VecStart_x)):
            #     ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'b')

            # VecStart_x = [-50, -50, -50, -50]
            # VecStart_y = [-5, -5, 5, 5]
            # VecStart_z = [0, 0, 0, 10]
            # VecEnd_x = [-50, -50, -50, -50]
            # VecEnd_y = [-5, 5, 5, -5]
            # VecEnd_z  =[10, 0, 10, 10]
            # for i in range(0, len(VecStart_x)):
            #     ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'b')

            # # t = 100
            # # t = (medida_z - random_plane_point[2][0]) / normal_Vec[2][0]
            # t = (-10 - random_plane_point[2][0]) / normal_Vec[2][0]
            # VecStart_x = [random_plane_point[0][0]]
            # VecStart_y = [random_plane_point[1][0]]
            # VecStart_z = [random_plane_point[2][0]]
            # VecEnd_x = [random_plane_point[0][0] + normal_Vec[0][0] * t]
            # VecEnd_y = [random_plane_point[1][0] + normal_Vec[1][0] * t]
            # VecEnd_z  =[random_plane_point[2][0] + normal_Vec[2][0] * t]

            # ax.scatter(VecStart_x[0], VecStart_y[0], VecStart_z[0], color = 'k')

            # for i in range(1):
            #     ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'k')
            
            # ax.scatter(0,0,0, color = 'k')
            # ax.scatter(x_1[0], y_1[0], medida_z, color = 'r')
            # ax.scatter(x_2[0], y_2[0], 0, color = 'r')
            # ax.scatter(medida_x, y_3[0], z_3[0], color = 'g')
            # ax.scatter(-medida_x, y_4[0], z_4[0], color = 'g')
            # ax.scatter(x_5[0], medida_y, z_5[0], color = 'purple')
            # # ax.scatter(random_plane_point[0], random_plane_point[1], random_plane_point[2], 'k')
            # # ax.plot([random_plane_point[0], 0], [random_plane_point[1],0], [random_plane_point[2],0])
            # # ax.scatter(x_2, y_2,0, 'k')

            # plt.show()


            list_z = [z_3, z_4, z_5, z_6]
            # print('Valores de t: ', t_1[0], t_2[0], t_3[0], t_4[0], t_5[0], t_6[0])
            # print('Valores de x: ', x_1[0], x_2[0], x_5[0], x_6[0])
            # print('Valores de y: ', y_1[0], y_2[0], y_3[0], y_4[0])
            # print('Valores de z: ', z_3[0], z_4[0], z_5[0], z_6[0])
            
            if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
                flag_cara_1 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
                flag_cara_2 = True
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
                flag_cara_3 = True
                # print('Bandera 3: ', flag_cara_3) 

            if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
                flag_cara_4 = True
                # print('Bandera 4: ', flag_cara_4)

            if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
                flag_cara_5 = True
                # print('Bandera 5: ', flag_cara_5) 

            if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
                flag_cara_6 = True
                # print('Bandera 6: ', flag_cara_6)

            list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
            
            Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

            if Delta_L > 0:
                # if Delta_L > 0 and Delta_L < 2.1:
                    list_delta_L.append(Delta_L)
                    list_random_th.append(np.degrees(Random_th[0]))
                    list_random_phi.append(np.degrees(Random_phi))

                    n_muons_in_CCD = n_muons_in_CCD + muon

                # else:
                #     n_negative_long = n_negative_long + 1
                #     continue

            else:
                    continue
            
    return list_delta_L, n_muons_in_CCD, n_negative_long, list_random_th, list_random_phi

def deffuse_model(number_thet,Theta, Theta_true, Phi, Radio, number_points_per_angle, long_a, long_b, 
                        medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z, dict_diffuse_model):
    
    list_P_vector = []
    list_random_point = []
    list_delta_L = []
    n_muons_in_CCD = 0
    n_negative_long = 0

    for i in np.arange(0,number_thet):
        list_random_energy_per_plane = []
        list_points_per_plane = []

        Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true
        Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi

        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # print(type(Vec[0]))
        Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        # norma = np.sqrt(Point[0] ** 2 + Point[1] ** 2 + Point[2] ** 2)
        # Point = (Vec[0], Vec[1], Vec[2])  ## Genera un punto sobre la esfera.
        # print('Vector sobre la esfera: ', Point)
        # print('Norma del Vector:', norma)
        # Points.append(Point)

        normal_Vec = (-1 * Vec[0] / Norma, -1 * Vec[1] / Norma, -1 * Vec[2] / Norma)     ## Es un vector normal unitario apuntando 
                                                                                            ##  hacia el centro de coordenadas
        # normal_Norma_Vec = norma_vec(normal_Vec)
        # print('Norma del vector anti-normal a la esfera:', normal_Norma_Vec)
        # print(len(normal_Vec))
        # Vectors.append(normal_Vec)

        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), np.cos(Random_th) * np.sin(Random_phi), np.sin(Random_th)]
        vec_phi = [-np.sin(Random_phi), np.cos(Random_phi), 0]
        # print('Vector Unitario Theta: ', vec_thet)
        # print('Vector Unitario Theta: ', vec_phi)

        for i in np.arange(0,number_points_per_angle):
            flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False, False, False, False 
            # list_random_th.append(Random_th[0])    ## Lo anexa en una lista
            # list_random_phi.append(Random_phi)


            random_a = rand.choice(long_a)  ## Selecciona un valor uniforme para el parámetro a
            random_b = rand.choice(long_b)  ##      ''      ''      ''      ''          ''    b

            # list_random_a.append(random_a)
            # list_random_b.append(random_b)

            P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                        random_a * vec_thet[1] + random_b * vec_phi[1], 
                        random_a * vec_thet[2] + random_b * vec_phi[2]]
            
            list_P_vector.append(P_vector)

            random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]
            # random_plane_point = [-1 * (Point[0] + P_vector[0]), -1 * (Point[1] + P_vector[1]), -1 * (Point[2] + P_vector[2])]

            # print(random_plane_point)
            list_random_point.append(random_plane_point)

            #### Intersecciones con cada cara   ####


            #### Cara Superior ###
            t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
            x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
            y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

            #### Cara Inferior ###
            t_2 = (0 - random_plane_point[2]) / normal_Vec[2] 
            x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
            y_2 = random_plane_point[1] + normal_Vec[1] * t_2

            ### Caras en X ###
            ### Cara 3 ###
            t_3 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
            z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
            y_3 = random_plane_point[1] + normal_Vec[1] * t_3

            ### Cara 4 ###
            t_4 = (medida_x - random_plane_point[0]) / normal_Vec[0]
            z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
            y_4 = random_plane_point[1] + normal_Vec[1] * t_4

            #### Caras en Y ###
            ### Cara 3 ###
            t_5 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
            z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
            x_5 = random_plane_point[0] + normal_Vec[0] * t_5

            ### Cara 4 ###
            t_6 = (medida_y - random_plane_point[1]) / normal_Vec[1]
            z_6 = random_plane_point[2] + normal_Vec[2] * t_4 
            x_6 = random_plane_point[0] + normal_Vec[0] * t_4

            list_z = [z_3, z_4, z_5, z_6]
            
            if np.around(x_1[0], 4) in mapeo_x and np.around(y_1[0], 4) in mapeo_y:
                flag_cara_1 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_2[0], 4) in mapeo_x and np.round(y_2[0], 4) in mapeo_y:
                flag_cara_2 = True
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3[0], 4) in mapeo_y and np.around(z_3[0], 4) in mapeo_z:
                flag_cara_3 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(y_4[0], 4) in mapeo_y and np.round(z_4[0], 4) in mapeo_z:
                flag_cara_4 = True
                # print('Bandera 2: ', flag_cara_2)

            if np.around(x_5[0], 4) in mapeo_x and np.around(z_5[0], 4) in mapeo_z:
                flag_cara_5 = True
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_6[0], 4) in mapeo_x and np.round(z_6[0], 4) in mapeo_z:
                flag_cara_6 = True
                # print('Bandera 2: ', flag_cara_2)

            list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
            
            Delta_L, muon = intersection_CCD(list_flags, list_z, medida_z, Random_th)

            if Delta_L != 0:
                if Delta_L > 0 and Delta_L < 2.1:
                    list_delta_L.append(Delta_L) 

                    ### -------------------- Modelo de Difusión ----------------- ##

                    # CCD_deep = 725 # micras
                    # Alpha = dict_diffuse_model['Alpha']
                    # Beta = dict_diffuse_model['Beta']

                    # z_deeps = np.arange()





                    n_muons_in_CCD = n_muons_in_CCD + muon

                else:
                    n_negative_long = n_negative_long + 1
                    continue

            else:
                    continue
            
    return list_delta_L, n_muons_in_CCD, n_negative_long

### ========================== Funciones de la pixelización ============================= ###
def diffution_curve(z, alpha, beta):
    return np.sqrt(np.abs(alpha * np.log(1 - (beta * z))))

def pixelizacion(alpha, beta, h, edep, theta, L):
    ratio_micra_cm = 10000 # micra / cm
    pixel_size = 15 # micras/px

    ### ------------------ Parámetros extra para dibujar la traza ----------------- ###
    CCD_thickness = 725 # micras

    CCD_thickness_px = CCD_thickness / pixel_size # px

    ### --------------------- Variables -------------------------- ###
    muon_height = h * ratio_micra_cm # micras
    muon_height_px = muon_height / pixel_size

    ratio_eV_electron = 3.7 # eV/e-

    Edep = edep # KeV
    Edep_eV = Edep * 10**(3) # eV
    Edep_elec = Edep_eV / ratio_eV_electron # e-

    Thet = np.degrees(theta) # Degrees
    Delta_L = L * ratio_micra_cm # micras
    ### ----------------------------------------------------------- ###


    delta_XY = np.sqrt(Delta_L**2 - muon_height**2 ) # micras
    # delta_XY = Delta_L * np.sin(Thet) # micras

    delta_XY_px = delta_XY / pixel_size # px
    list_XY_deep = np.linspace(0, muon_height_px, int(np.around(delta_XY_px, 0)))

    z_deeps = np.arange(7.5, CCD_thickness, 7.5) # micras

    sigma_values = diffution_curve(list_XY_deep * pixel_size, alpha, beta) # micras
    max_width = int(np.around(sigma_values[-1] * 4 , 0))


    ### Distribución de la carga por linea ###
    charge_per_line = Edep_elec / delta_XY_px



    return

#### ======================== FUNCIONES PARA LA BARRA CENTELLADORA ============================== #### 
def LandV_Barr(lx, lpar):
    Delta = lx[0]	# Energy loss in absorber
    L = lpar[0]		# Thickness of absorber (Distance crossed by the particle)
    p = lpar[1]		# Momentum (in MeV/c)

    K = 0.1535	# K coefficient = 2*pi*N*r^2*m*c^2 (in MeV g^-1 cm^2)

    z = -1;		# Charge number of incident particle
    ZA = 0.54141 # Atomic number over Atomic mass of absorber (for PVT in this case)
    c = TMath.C()	# Speed of light
    me = 0.510998928	#  Electron mass in MeV/c^2
    M = 105.65839	# Muon mass in MeV/c^2
    I = 64.7/1000000; # Mean excitation energy in MeV (for PVT)

    bg = p/M
    beta = bg/np.sqrt(1+bg**2)	# Beta factor
    gamma = 1/np.sqrt(1-beta**2)	# Gamma factor
    pi = TMath.Pi()
    rho = 1.032 # Density of material in gr/cm^3(for PVT)

    a = 0.1610		# Parameters (taken from PDG for PVT)
    k = 3.2393			
    X0 = 0.1464
    X1 = 2.4855
    C = 3.1997
    d0 = 0.0	

    X = np.log10(bg)
    # d;	// Variable for the Density effect
    
    if X>=X1:
        d = 2 * np.log(10.0) * X - C
    elif X0<=X and X<X1:
        d = 2 * np.log(10.0) * X - C + a * ((X1-X)**k)
    elif X<X0:
        d = d0 * 10**(2*(X-X0))

    WM = 2 * me * ((bg)**2)/(1+(2*me*gamma/M)+(me/M)**2)   #Maximum energy tranfer (from PDG)
    # WM = (2 * me * (beta*gamma)**2)/(1+(2*me/M)*np.sqrt(1 + (beta*gamma)**2) + (me/M)**2)   #Maximum energy tranfer

    loge = np.log((1-beta**2)*(I**2)/(2*me*(beta**2))) + (beta**2) # log epsilon variable (with electron mass)
    # loge = np.log((1-beta**2)*(I**2)/(2*M*(beta**2))) + (beta**2) # log epsilon variable (with muon mass)

    EC = 0.577	# Euler's constant

    xi = (K)*rho*ZA*L*(z/beta)**2		# Xi variable 
    DeltaAv = xi * (np.log(2*me*(gamma**2)*(beta**2)*WM/(I**2))-(2*beta**2)-(d)) # Mean energy loss (Bethe-Bloch of PDG)
    # DeltaAv = xi * (np.log(2*me*(gamma**2)*(p**2)*WM/(I**2))-(2*beta**2)-(d)) # Mean energy loss (Bethe-Bloch)

    Lambda = (Delta-xi*(np.log(xi)-loge+1-EC))/xi # Lambda parameter
    # Lambda = (Delta - DeltaAv)/xi - beta**2 - np.log(xi/I) - 1 + EC # Lambda parameter (PDG?)

    # Deltamp = xi*(np.log(xi/np.exp(loge))+0.198-d)		# Most probable energy loss (Leo)
    Deltamp = xi*(np.log(xi) - loge + 0.198 - d)		# Most probable energy loss (Leo)
    lambdamp = (Deltamp-xi*(np.log(xi)-loge+1-EC))/xi 

    # Deltamp = xi * (np.log((2 * me * beta**2 * gamma**2)/I) + np.log(xi/I) + 0.2 - beta**2 - d)		# Most probable energy loss from PDG
    # lambdamp = (Deltamp-xi*(np.log(xi)-loge+1-EC))/xi

    kappa = xi/WM		# Kappa ratio
    beta2 = beta**2
    sigma2 = (xi**2)*(1-beta2/2)/kappa		# Standard deviation for relativistic particles

    if kappa<=0.01:
        # phi = TMath.Landau(Lambda, Deltamp, 1.0)
        phi = TMath.Landau(Lambda, lambdamp, 1.0)
        return phi/xi

    # elif kappa>0.01 and kappa<10:
    elif 0.01<kappa and kappa<10:
        # vav = TMath.Vavilov(Delta-Deltamp, kappa, beta2)
        vav = TMath.Vavilov(Lambda - lambdamp, kappa, beta2)
        return vav

    else:
        gauss = TMath.Gaus(Delta, DeltaAv, np.sqrt(sigma2))
        return gauss

def random_LV_Barr(s, p):
    gRandom.SetSeed(0) ## Cambia la semilla aleatoria para el GetRandom
    f = TF1("f", LandV_Barr, 0, 100, 2)
    # f = TF1("", LandV, 0, 10, 2)
    # f.SetParameter(0, np.double(s))
    # f.SetParameter(1, np.double(p))
    f.SetParameters(s, p)

    Edep = f.GetRandom()

    K = 0.1535	# K coefficient = 2*pi*N*r^2*m*c^2 (in MeV g^-1 cm^2)

    z = -1;		# Charge number of incident particle
    ZA = 0.54141 # Atomic number over Atomic mass of absorber (for PVT in this case)
    c = TMath.C()	# Speed of light
    me = 0.510998928	#  Electron mass in MeV/c^2
    M = 105.65839	# Muon mass in MeV/c^2
    I = 64.7/1000000; # Mean excitation energy in MeV (for PVT)

    bg = p/M
    beta = bg/np.sqrt(1+bg**2)	# Beta factor
    gamma = 1/np.sqrt(1-beta**2)	# Gamma factor
    pi = TMath.Pi()
    rho = 1.032 # Density of material in gr/cm^3(for PVT)

    xi = (K)*rho*ZA* s *(z/beta)**2		# Xi variable
    WM = 2 * me * ((bg)**2)/(1+(2*me*gamma/M)+(me/M)**2)   #Maximum energy tranfer (Bryan from PDG)

    kappa = xi/WM		# Kappa ratio


    return Edep * 1000, kappa    #### Energía en KeV

def muon_generator_BARRA(number_thet, Radio, medida_x, medida_y, medida_z, mapeo_x, mapeo_y, mapeo_z, half_plane_size):
    list_thet_in_CCD = []
    list_phi_in_CCD = []
    list_energy_pri_in_CCD = []
    list_nmuons = []
    list_delta_L = []
    list_energy_Landau = []
    list_kappa = []

    nmuons_in_CCD = 0
    n_negative_long = 0
    muon_in_bucle = 0

    max_long = np.sqrt((medida_x * 2)**2 + (medida_y * 2)**2 + (medida_z * 2)**2)

    for i in np.arange(0,number_thet):
        # n_muon = i
        # Random_th = rand.choices(Theta, Theta_true) ## Escoje un ángulo segun la distribución de Theta_true en radianes
        # Random_phi = rand.choice(Phi)   ## Lo mismo pero con phi en radianes

        Random_th = random_thet() ## Escoje un ángulo segun la distribución de Theta_true en radianes
        Random_phi = rand.random() * (2 * np.pi) ## Lo mismo pero con phi en radianes
        Random_energy = random_SD(Random_th) ## Escoje una energía segun la distribución de Smith-Duller en MeV

        # print(np.cos(Random_th))
        # print(Random_th)
        # print(Random_phi)
        # print(Random_energy)

        ### ==================== Momento del muon ======================= ###
        m_mu = 105.7    ## MeV/c^2
        En_tot = Random_energy + m_mu  
        momentum = np.sqrt(En_tot**2 - m_mu**2)
        ### ============================================================= ###

        Vec = coord_cartesian(Random_th, Random_phi)
        Norma = norma_vec(Vec)
        # Point = [Radio * Vec[0], Radio * Vec[1], Radio * Vec[2]]  ## Genera un punto sobre la esfera.
        Point = [Radio * np.sin(Random_th) * np.cos(Random_phi), 
                 Radio * np.sin(Random_th) * np.sin(Random_phi), 
                 Radio * np.cos(Random_th)]

        normal_Vec =  [-1 * np.sin(Random_th) * np.cos(Random_phi), 
                       -1 * np.sin(Random_th) * np.sin(Random_phi), 
                       -1 * np.cos(Random_th)]     ## Es un vector normal unitario apuntando hacia el centro de coordenadas

        ### ======= Vectores unitarios ====== ###
        vec_thet = [np.cos(Random_th) * np.cos(Random_phi), 
                    np.cos(Random_th) * np.sin(Random_phi),
                    -np.sin(Random_th)]
        
        vec_phi = [-np.sin(Random_phi), 
                   np.cos(Random_phi), 
                   0]
        ### ================================= ###

        random_a = -half_plane_size + rand.random()*2 * half_plane_size ## Selecciona un valor uniforme para el parámetro a
        random_b = -half_plane_size + rand.random()*2 * half_plane_size ##      ''      ''      ''      ''          ''    b

        P_vector = [random_a * vec_thet[0] + random_b * vec_phi[0], 
                    random_a * vec_thet[1] + random_b * vec_phi[1], 
                    random_a * vec_thet[2] + random_b * vec_phi[2]]

        random_plane_point = [Point[0] + P_vector[0], Point[1] + P_vector[1], Point[2] + P_vector[2]]

        #### Intersecciones con cada cara   ####
        flag_cara_1, flag_cara_2, flag_cara_3 = False, False, False
        flag_cara_4, flag_cara_5, flag_cara_6 = False, False, False

        #### Cara Superior 1 ###
        t_1 = (medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_1 = random_plane_point[0] + normal_Vec[0] * t_1 
        y_1 = random_plane_point[1] + normal_Vec[1] * t_1 

        #### Cara Inferior 2 ###
        t_2 = (-medida_z - random_plane_point[2]) / normal_Vec[2] 
        x_2 = random_plane_point[0] + normal_Vec[0] * t_2 
        y_2 = random_plane_point[1] + normal_Vec[1] * t_2

        ### Caras en X ###
        ### Cara 3 ###
        t_3 = (medida_x - random_plane_point[0]) / normal_Vec[0]
        z_3 = random_plane_point[2] + normal_Vec[2] * t_3 
        y_3 = random_plane_point[1] + normal_Vec[1] * t_3

        ### Cara 4 ###
        t_4 = (-medida_x - random_plane_point[0]) / normal_Vec[0]
        z_4 = random_plane_point[2] + normal_Vec[2] * t_4 
        y_4 = random_plane_point[1] + normal_Vec[1] * t_4

        #### Caras en Y ###
        ### Cara 5 ###
        t_5 = (medida_y - random_plane_point[1]) / normal_Vec[1]
        z_5 = random_plane_point[2] + normal_Vec[2] * t_5 
        x_5 = random_plane_point[0] + normal_Vec[0] * t_5

        ### Cara 6 ###
        t_6 = (-medida_y - random_plane_point[1]) / normal_Vec[1]
        z_6 = random_plane_point[2] + normal_Vec[2] * t_6
        x_6 = random_plane_point[0] + normal_Vec[0] * t_6

        list_z = [z_3, z_4, z_5, z_6]
        n_flags = 0
        limit_around = 3
        
        flag_faces = True
        while flag_faces:
            if np.around(x_1, limit_around) in mapeo_x and np.around(y_1, limit_around) in mapeo_y:
                flag_cara_1 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_2, limit_around) in mapeo_x and np.round(y_2, limit_around) in mapeo_y:
                flag_cara_2 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)

            if np.around(y_3, limit_around) in mapeo_y and np.around(z_3, limit_around) in mapeo_z:
                flag_cara_3 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(y_4, limit_around) in mapeo_y and np.round(z_4, limit_around) in mapeo_z:
                flag_cara_4 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)

            if np.around(x_5, limit_around) in mapeo_x and np.around(z_5, limit_around) in mapeo_z:
                flag_cara_5 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 1: ', flag_cara_1) 

            if np.round(x_6, limit_around) in mapeo_x and np.round(z_6, limit_around) in mapeo_z:
                flag_cara_6 = True
                n_flags = n_flags + 1
                if n_flags == 2:
                    flag_faces = False
                # print('Bandera 2: ', flag_cara_2)
            else:
                flag_faces = False

        list_flags = [flag_cara_1, flag_cara_2, flag_cara_3, flag_cara_4, flag_cara_5, flag_cara_6]
        
        Delta_L, _ = intersection_CCD(list_flags, list_z, medida_z, Random_th)

        if Delta_L < 0:
            print(Delta_L)

        # Fin = datetime.datetime.now()
        # print('Tiempo de cálculo para Delta L: ', Fin-In)

        if 0 < Delta_L and Delta_L < max_long:
            # print(list_flags)
            # print(list_z)
            nmuons_in_CCD =  nmuons_in_CCD + 1
            # print("Momentum: ", momentum, ", Distance: ", Delta_L)

            # try:
            Random_energy_Landau, kappa = random_LV_Barr(s=Delta_L, p = momentum)
            # print('Energy Landau: ', Random_energy_Landau)
            print('Kappa: ', kappa)

            # except: 
            #     print('Hubo en error con el cálculo de Landau')
            #     print(Random_energy_Landau, kappa)
            #     # break
            #     continue


            # list_nmuons.append(n_muon)
            list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)
            list_kappa.append(kappa)

            # list_rand_thet.append(Random_th[0])
            # list_rand_phi.append(Random_phi)
            
            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '\r')

            # else:
            #     n_negative_long = n_negative_long + 1
            #     continue

        else:
            # print('DeltaL negativo')
            Random_energy_Landau = 0 # En KeV
            Delta_L = 0
            kappa = -10

            # list_nmuons.append(n_muon)
            list_energy_Landau.append(Random_energy_Landau)
            list_thet_in_CCD.append(Random_th)
            list_phi_in_CCD.append(Random_phi)
            list_energy_pri_in_CCD.append(Random_energy)
            list_delta_L.append(Delta_L)
            list_kappa.append(kappa)

            muon_in_bucle += 1

            print('Muon simulado ' + str(muon_in_bucle) + '/' + str(number_thet), end = '\r')

    dict_muons =  {'NMuon': list_nmuons, 'Theta(Rad)': list_thet_in_CCD, 'Phi(Rad)' : list_phi_in_CCD, 'Energy-SD(MeV)' : list_energy_pri_in_CCD, 
                        'Delta_L(cm)' : list_delta_L, 'Energy_Landau(KeV)' : list_energy_Landau, 'kappa' : list_kappa} 

    return dict_muons, nmuons_in_CCD 

def dimension_x_barr(long_x):
    step = 0.001

    list_long_x = [-long_x]

    while long_x:
        x = np.round(list_long_x[-1] + step, 3)
        list_long_x.append(x)

        if list_long_x[-1] == long_x:
            break

    return list_long_x

def dimension_y_barr(long_y):
    step = 0.001

    list_long_y = [-long_y]

    while long_y:
        y = np.round(list_long_y[-1] + step, 3)
        list_long_y.append(y)

        if list_long_y[-1] == long_y:
            break

    return list_long_y

def dimension_z_barr(long_z):
    step = 0.001

    list_long_z = [-long_z]
    # list_long_z = [0]

    while long_z:
        z = np.round(list_long_z[-1] + step, 3)
        list_long_z.append(z)

        if list_long_z[-1] == long_z:
            break

    return list_long_z

#### ============================================================================================ ####

