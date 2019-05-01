# -*- coding: utf-8 -*-
"""
@Title:     Physics 660 Project One
@Author     Andrew Crossman
@Date       Feb. 12th, 2018
"""
#Imports
import numpy as np
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Static Parameters
initial_a = 200
initial_b = 5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Part One Equation Definition
def part1(tau_a, tau_b, dt):
    #initialize time steps and data arrays
    time_steps = int(10/dt)
    Na = np.zeros(time_steps)
    exact_Na = np.zeros(time_steps)
    Nb = np.zeros(time_steps)
    exact_Nb = np.zeros(time_steps)
    #allocate initial variables
    Na[0] = initial_a
    exact_Na[0] = initial_a
    Nb[0] = initial_b
    exact_Nb[0] = initial_b
    time = np.zeros(time_steps)
    #Determine analytical and numerical calculations
    for i in list(range(time_steps-1)):
        time[i+1] = time[i]+dt
        #N_a(t) equations
        Na[i+1] = Na[i] - (Na[i]/tau_a) * dt
        exact_Na[i+1] = initial_a*np.exp(-time[i+1]/tau_a)
        #N_b(t) equations
        Nb[i+1] = Nb[i] - (Nb[i]/tau_b) * dt + (Na[i]/tau_a) * dt
        #use line below when tau_a =/= tau_b
        if tau_a!=tau_b:
            exact_Nb[i+1] = ((initial_a * tau_b * np.exp(-time[i+1]/tau_a))/(tau_a - tau_b)) + ((5 - (200*tau_b/(tau_a - tau_b))) * np.exp(-time[i+1]/tau_b))
        #use line below for when tau_a = tau_b
        else:
            exact_Nb[i+1] = (initial_a * time[i+1] *np.exp(-time[i+1]/tau_a))/(tau_a) + (initial_b*np.exp(-time[i+1]/tau_a))   
    #Plot numeric and analytic solutions
    f,ax = plt.subplots()
    ax.scatter(time,Na,marker='*',color='k',label='Na(t)')
    ax.plot(time,exact_Na,color='r',label='exact Na(t)')
    ax.scatter(time,Nb,marker='p',color='c',label='Nb(t)')
    ax.plot(time,exact_Nb,color='g',label='exact Nb(t)')
    if tau_a==tau_b:
        ax.set_title("Number of Radioactive Nuclei: "+r'$\tau_a = \tau_b$' ,style='italic')        
        ax.annotate('Equal Populations', xy=(195*tau_a/200, 73), xytext=(3, 125),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
    elif tau_a>tau_b:
        ax.set_title("Number of Radioactive Nuclei: "+r'$\tau_a > \tau_b$' ,style='italic')
    else:
        ax.set_title("Number of Radioactive Nuclei: "+r'$\tau_a < \tau_b$' ,style='italic')
        xcord =np.log((200-200*tau_b/(tau_a-tau_b))/(5-200*tau_b/(tau_a-tau_b)))*(1)/((1/tau_a)-(1/tau_b))
        ycord =200*np.exp(-xcord/tau_a)
        print(xcord)
        ax.annotate('Equal Populations', xy=(xcord, ycord), xytext=(3, 125),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
    if tau_a!=tau_b:
        ax.set_xlabel("Time ("+r'$t/\tau_a$'+")",style='italic')
    else:
        ax.set_xlabel("Time ("+r'$t/\tau$'+")",style='italic')
    ax.set_ylabel("Number of Nuclei",style='italic')
    ax.legend(loc='upper right')
    ax.annotate(r'$\tau_a=$'+str(tau_a)+'\n'+r'$\tau_b=$'+str(tau_b),(2,180))
    if tau_a==tau_b:
        f.savefig("Radioactive_DecayAequalB.png",dpi=600)
    elif tau_a>tau_b:
        f.savefig("Radioactive_DecayAgreaterthanB.png",dpi=600)
    else:
        f.savefig("Radioactive_DecayAlessthanB.png",dpi=600)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Part Two Equation Defintion
def part2(tau, dt):
    #initialize time steps and data arrays
    time_steps = int(10/dt)
    #time_steps = int(10*tau/dt)
    Na = np.zeros(time_steps)
    exact_Na = np.zeros(time_steps)
    Nb = np.zeros(time_steps)
    exact_Nb = np.zeros(time_steps)
    #allocate initial values
    Na[0] = 150
    exact_Na[0] = 150
    Nb[0] = 0
    exact_Nb[0] = 0
    time = np.zeros(time_steps)
    #Determine analytical and numerical calculations
    for i in list(range(time_steps-1)):
        time[i+1] = time[i]+dt
        
        Na[i+1] = 50*(1+2*np.exp((-3*time[i+1])/(2*tau))) + ((150/tau)*np.exp((-3*time[i+1])/(2*tau)))*dt
        exact_Na[i+1] = 50 + 100*np.exp((-3*time[i+1])/(2*tau))
        
        Nb[i+1] = 100*(1-np.exp((-3*time[i+1])/(2*tau))) - ((150/tau)*np.exp((-3*time[i+1])/(2*tau)))*dt
        exact_Nb[i+1] = 100 - 100*np.exp((-3*time[i+1])/(2*tau))
        
    f,ax = plt.subplots()
    ax.scatter(time,Na,marker='*',color='k',label='Na(t)')
    ax.plot(time,exact_Na,color='r',label='exact Na(t)')
    ax.scatter(time,Nb,marker='p',color='c',label='Nb(t)')
    ax.plot(time,exact_Nb,color='g',label='exact Nb(t)')
    ax.annotate('Equal Populations', xy=(2*tau*np.log(4)/3, 75), xytext=(3, 140),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
    ax.set_title("Number of Radioactive Nuclei: "+r'$\tau = $'+ str(tau) ,style='italic')        
    ax.set_xlabel("Time ("+r'$t/\tau$'+")",style='italic')
    ax.set_ylabel("Number of Nuclei",style='italic')
    ax.legend(loc='upper right')
    f.savefig("Radioactive_DecayPart2.png",dpi=600)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
#Main Code
#Calls equation part1 for the three cases of tau_a and tau_b
part1(1,1,.25)
part1(1,5,.25)
part1(1,.2,.25)
#Calls eqaution for part2
part2(1,.25)