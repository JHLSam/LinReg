import numpy as np
import matplotlib.pyplot as plt
import math
#from mpl_toolkits.mplot3d import Axes3D
plt.rcParams.update({'errorbar.capsize':2}) #error with rcparams in current version of API, errorbar caps will not show up without this update
#!/usr/bin/env python
# -*- coding: utf-8 -*-

y = (4.1268*10**-16,5.05*10**-16,1.769633*10**-15,3.36563*10**-15,5.74063*10**-15,9.54959*10**-15,1.41631*10**-14,2.21054*10**-14,2.40622*10**-14,3.26904*10**-14,3.5879*10**-14,4.78998*10**-14) #create array(list) w/ range 0 to 10(exclusive) spacing:0.1
x = (996,9.97*10**3,99.6*10**3,215.2*10**3,386*10**3,660*10**3,0.992*10**6,1.486*10**6,1.744*10**6,2.244*10**6,2.557*10**6,3.38*10**6)
less2df = len(x)-2 #n-2 degrees of freedom
reciprocal_df = 1/len(x) #degrees of freedom reciprocal
reciprocal_dfless2 = 1/less2df #n-2 degrees of freedom reciprocal
mx = np.mean(x)
mx_squared = np.mean(x)**2
my = np.mean(y)
denom = x - mx
c1 = np.sum(denom*(y-my))/np.sum(denom ** 2) #gradient
c0 = my - c1 * mx #y-intercept

def scalarprod2(k):

    x = [996,9.97*10**3,99.6*10**3,215.2*10**3,386*10**3,660*10**3,0.992*10**6,1.486*10**6,1.744*10**6,2.244*10**6,2.557*10**6,3.38*10**6]
    z = []
    for i in x:
        z.append(i*c1 * k)
    return z

z1 = scalarprod2(1)
z2 = np.array(z1)
d1 = y - z2 - c0
slope_estimate = math.sqrt(np.sum(d1**2)/np.sum(denom**2)*(reciprocal_df)) #gradient uncertainty
print("slope:",c1)
print("slope estimate:", slope_estimate)
t1 = reciprocal_df + np.mean(x)**2/np.sum(denom**2)
t2 = np.sum(d1**2)*reciprocal_dfless2
intercept_estimate = math.sqrt(t1*t2) #y-intercept uncertainty  
print("y-intercept:",c0)
print("intercept_estimate:", intercept_estimate)
y2 = [0,4*10**6] #x value interval
x2 = [c0 + c1*0, c0 + c1*4*10**6] #y value interval 
z2 = [0,100]

def scalprod(x):

    y = [4.1268*10**-16,5.05*10**-16,1.769633*10**-15,3.36563*10**-15,5.74063*10**-15,9.54959*10**-15,1.41631*10**-14,2.21054*10**-14,2.40622*10**-14,3.26904*10**-14,3.5879*10**-14,4.78998*10**-14]
    y3 = []
    for i in y:
        y3.append(x*i)
    return y3

y4 = scalprod(0.05) #5% y-value uncertainty
y5 = np.array(y4) #transform list in scalprod to array
my_dpi = 96
fig1=plt.figure(figsize=(800/my_dpi, 600/my_dpi), dpi=my_dpi)
plt.plot(y2,x2, color='black', linewidth=3)
plt.suptitle("⟨V⟩"+'$^{2}$/G'+'$_{0}$'+'$^{2}$'+'f'+'$_{0}$'+'π/2'+"(1/hz) vs Resistance\nName:\nVersion:Thermal Noise\n")
#plt.axis([0,2*10**7,0,2])
plt.errorbar(x,y,xerr=None, yerr=0.08*10**-14,fmt='o',capsize=2)
plt.text(2500000,2*10**-14,'m±∆m:(1.4*E-20.1±1.1*E-22)\nc±∆c:3.6*E-16±2*E-16',fontsize=10)
#plt.text(0.2,80,'m+/-(del)m:'+ slope_estimate + c+/-(del)c:0.7+/-0.2',fontsize=10)
plt.xlabel("Resistance(Ω)")
plt.ylabel("⟨V⟩"+'$^{2}$/G'+'$_{0}$'+'$^{2}$'+'f'+'$_{0}$'+'π/2'+"(1/hz)")
plt.savefig("lin_regression(b).png")
plt.show(fig1)
