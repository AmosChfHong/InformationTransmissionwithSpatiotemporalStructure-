# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 16:39:06 2015

@author: Hong
"""

#%% General Descriptions
from brian2 import *
import SpikeAnalysis as mi
import numpy as np
import copy

# Control Variables
defaultclock.dt = 10*us
SimulationTime = 100*ms
PatternNum = 4#64
RepTimes = 20#64
InputNum = 400
TestNum = 40

#%% Generate Input
InputPatterns = []
start_scope()

t0 = 10*ms
si = 3*ms
am = 0.55
for i in range(PatternNum):
    print("Generating %dth pattern in %d patterns"%((i+1),PatternNum))
    InputNeuron =  NeuronGroup(InputNum,model="U=am*exp(-(t-t0)**2/(2*si**2))/((2*pi)**0.5*si/ms):1", 
                      threshold='rand()<U*dt/ms')                   
    InputMon= SpikeMonitor(InputNeuron)
    run(SimulationTime)
    InputIndex = []
    InputTimes = []     
    InputTrain = InputMon.values('t')
    for i in range(InputNum):
        for T in InputTrain[i]:
            InputIndex.append(i)
            InputTimes.append(T+t0-si)    
    InputPatterns.append([InputIndex,InputTimes])         


    



#%% Neuron Construction
start_scope()
#-------------Parameters--------------#
#Neuron Parameters
#Cm = 400*pF
Veq = -55*mV
Vre = -65*mV
Vth = -45*mV
gL = 45*nS
Ee = 0*mV
Ei = -80*mV
tau_se = 2*ms
tau_si = 5*ms
ENoiseRat = 2*Hz
INoiseRat = 12.5*Hz
ENoiseNum = 18000*4
INoiseNum = 2000*4
#Neuron Group Parameters
ExcNum = InputNum
LatNum = 1
LayerNum = 5



#%% Neuron Construction
start_scope()



#-----------------------Input Source----------------------------------#
Index = []
Times = []*ms
InputNeuron = SpikeGeneratorGroup(InputNum,Index,Times)

#-----------------------Bulid Network Neurons--------------------------------#
#Differential equation for neuron dynamics
NeuDyn = '''dv/dt = (gL*(Veq-v) + gExc*(Ee-v) + gInh*(Ei-v))/Cm: volt  (unless refractory)
            dgExc/dt = -gExc/tau_se: siemens
            dgInh/dt = -gInh/tau_si: siemens
            Cm:farad
            '''

#Delare NeuronGroups                  
ExcNeuron = NeuronGroup(ExcNum*LayerNum, model=NeuDyn,
                        method='euler',threshold='v>Vth',
                        reset='v=Vre', refractory = 2*ms)                                               
LatNeuron = NeuronGroup(LatNum*(LayerNum),model=NeuDyn,
                        method='euler',threshold='v>Vth',
                        reset='v=Vre', refractory = 2*ms)
                                                                      
#Delare Independent Noise
w_e = 0.40*nS
w_i = 0.25*nS

#P1 = PoissonInput(ExcNeuron,'gExc',N=ENoiseNum,rate=ENoiseRat, weight=w_e)
#P2 = PoissonInput(ExcNeuron,'gInh',N=INoiseNum,rate=INoiseRat, weight=w_i)   
 
                                                                
P1 = PoissonGroup(N=ENoiseNum,rates=ENoiseRat)
P2 = PoissonGroup(N=INoiseNum,rates=INoiseRat)

P1Link = Synapses(P1, ExcNeuron, model="w:siemens",
                     on_pre="gExc_post += w")
P2Link = Synapses(P2, ExcNeuron, model="w:siemens",
                     on_pre="gInh_post += w")
P3Link = Synapses(P1, LatNeuron, model="w:siemens",
                     on_pre="gExc_post += w")
P4Link = Synapses(P2, LatNeuron, model="w:siemens",
                     on_pre="gInh_post += w")                     
P1Link.connect(True,p=0.25)
P2Link.connect(True,p=0.25)
P3Link.connect(True,p=0.25)
P4Link.connect(True,p=0.25)

P1Link.w[:,:] = w_e
P2Link.w[:,:] = w_i
P3Link.w[:,:] = w_e
P4Link.w[:,:] = w_i


ExcNeuron.v[:] = "-(47)*mV"
LatNeuron.v =  "-(46 + rand())*mV"
ExcNeuron.gExc = 10*nS
ExcNeuron.gInh = 10*nS
ExcNeuron.Cm = 400*pF
LatNeuron.Cm = 300*pF  
#%% Network Connection 

#------------------Parameters-----------------------------#
#Connection probability
P_ee = 0.2
P_el = 0.2
P_le = 0.6

#Snaptic Weight
w_ee = 0.9*nS #4.5*nS #2.0 ff
w_el = 00.0*nS 

w_le = 0.0*nS#3.6*nS

sigma_w = 0.4
#-----------------Connections----------------------------#
#Construct Syanpses 
IELink = Synapses(InputNeuron, ExcNeuron, model="w:siemens",
                     on_pre="gExc_post += w") # tau_se is used to normalize the input so that all conduction is equal to w

ILLink = Synapses(InputNeuron, LatNeuron, model="w:siemens",
                     on_pre="gExc_post += w") 

EELink = Synapses(ExcNeuron,ExcNeuron,model="w:siemens",
                  on_pre="gExc_post += w")   
                  
ELLink = Synapses(ExcNeuron,LatNeuron,model="w:siemens",
                  on_pre="gExc_post += w")  

LELink = Synapses(LatNeuron,ExcNeuron,model="w:siemens",
                  on_pre="gInh_post += w") 
                  
LLLink = Synapses(LatNeuron,LatNeuron,model="w:siemens",
                  on_pre="gInh_post += w")   
                  
#Declare Connections
IELink.connect('j<ExcNum', p = P_ee)  
ILLink.connect('j<LatNum', p = P_el)   
LELink.connect(''' i>=0 and i<LatNum and \
                   j>=0 and j<ExcNum''', p=P_le)   
LLLink.connect(''' i>=0 and i<LatNum and \
                   j>=0 and j<LatNum''', p=P_le)
for Ind in range(LayerNum-1):
    EELink.connect(''' i>=Ind*ExcNum and \
                      i<(Ind+1)*ExcNum and \
                      j>=(Ind+1)*ExcNum and \
                      j<(Ind+2)*ExcNum''', p = P_ee)
    ELLink.connect(''' i>=Ind*ExcNum and \
                      i<(Ind+1)*ExcNum and \
                      j>=(Ind+1)*LatNum and \
                      j<(Ind+2)*LatNum''', p=P_el)
    LELink.connect(''' i>=(Ind+1)*LatNum and \
                      i<(Ind+2)*LatNum and \
                      j>=(Ind+1)*ExcNum and \
                      j<(Ind+2)*ExcNum''', p=P_le)
    LLLink.connect(''' i>=(Ind+1)*LatNum and \
                      i<(Ind+2)*LatNum and \
                      j>=(Ind+1)*LatNum and \
                      j<(Ind+2)*LatNum''', p=P_le)
                      
    
#Set Weights
sigma = 2.0

LN = len(IELink.w)
IELink.w[:,:] = clip(np.random.lognormal(log(w_ee/nS)-0.5*sigma**2,sigma,LN)*nS,0,50*w_ee)
LN = len(ILLink.w)
ILLink.w[:,:] = clip(np.random.lognormal(log(w_el/nS)-0.5*sigma**2,sigma,LN)*nS,0,20*w_el)
LN = len(EELink.w)
EELink.w[:,:] = clip(np.random.lognormal(log(w_ee/nS)-0.5*sigma**2,sigma,LN)*nS,0,20*w_ee)
LN = len(ELLink.w)
ELLink.w[:,:] = clip(np.random.lognormal(log(w_el/nS)-0.5*sigma**2,sigma,LN)*nS,0,20*w_el)
LELink.w[:,:] = "clip(w_le*(1+sigma_w*randn()),0,5*w_le)"  
LLLink.w[:,:] = "clip(w_le*(1+sigma_w*randn()),0,5*w_le)"    
 
  
#hist(EELink.w,100),show()  
#Set Delays
delays = 6
#beta_d = 1
#LN = len(IELink.delay)
#IELink.delay[:,:] = (2 + delays*np.random.beta(beta_d,beta_d,LN))*ms
#LN = len(EELink.delay)
#EELink.delay[:,:] = (2 + delays*np.random.beta(beta_d,beta_d,LN))*ms
#LN = len(ILLink.delay)
#ILLink.delay[:,:] = (2 + delays*np.random.beta(beta_d,beta_d,LN))*ms
#LN = len(ELLink.delay)
#ELLink.delay[:,:] = (2 + delays*np.random.beta(beta_d,beta_d,LN))*ms
IELink.delay[:,:] = "(5+delays*rand())*ms"
EELink.delay[:,:] = "(5+delays*rand())*ms"
ILLink.delay[:,:] = "(5+delays*rand())*ms"
ELLink.delay[:,:] = "(2+delays*rand())*ms"
LELink.delay[:,:] = "(2*rand())*ms"
#%% Set Monitors
    
    
#Spike Monitors
SucSpk = SpikeMonitor(InputNeuron)    
ExcSpk = SpikeMonitor(ExcNeuron)   
LatSpk = SpikeMonitor(LatNeuron)    

#StateS Monitors     
ExcU = StateMonitor(ExcNeuron,('v','gExc','gInh'),record = True, dt=0.5*ms)      




store()

#%% Run Simulation
PopulationTimes = []
IndividualTimes = []
TrailPopulation = []
for i in range(LayerNum):
    PopulationTimes.append([])
    IndividualTimes.append([])
    TrailPopulation.append([])
    for j in range(ExcNum):
        IndividualTimes[i].append([])
               
print("runing the simulation...")                
Datas = []
ExcData = []
InpData = []

    
PatternData = []   
for ii in range(PatternNum):
    LayerTemp = []
    for i in range(LayerNum):       
        PopulationTemp = []
        for j in range(ExcNum):
              PopulationTemp.append([])
        LayerTemp.append(PopulationTemp)
            
        
    for jj in range(RepTimes):   
        restore()
        InputIndex = InputPatterns[ii][0]
        InputTimes = InputPatterns[ii][1]
        InputNeuron.set_spikes(array(InputIndex),array(InputTimes)*second)
        print("runing the %dth pattern %dth times"%((ii+1),(jj+1)))
        run(SimulationTime)
        #%% Data Collection 
        
        if(RepTimes > 1):  
            ExcTrain = ExcSpk.values('t')                    
            for i in range(LayerNum): 
                LayerT = []
                for j in range(ExcNum): 
                    for T in ExcTrain[j+ExcNum*i]:
                        LayerT.append(T/ms)
                Population = sort(LayerT)
                LP = len(Population)
                MeanTime = Population[LP/2]
                for j in range(TestNum): 
                    Times = []
                    for T in ExcTrain[j+ExcNum*i]:
                        Time = T/ms - MeanTime
                        if abs(Time)<10 :
                            Times.append(Time)
                    ExcData.append(Times) 
                    
            InpTrain = SucSpk.values('t')
            for Neuron in range(InputNum):
                SpikeTrain = sort(InputTrain[Neuron]/ms)
                Times = []
                for T in SpikeTrain:
                    Times.append(T)
                InpData.append(Times)
                
                
        ExcTrain = ExcSpk.values('t')                    
        for i in range(LayerNum):
            Temp = []
            for j in range(ExcNum):
                for T in ExcTrain[i*ExcNum + j]:
                    Temp.append(T/ms)
            Population = sort(Temp)
            LP = len(Population)
            if LP>0:
                MeanTime = Population[LP/2]
                for k in range(ExcNum):
                    for T in ExcTrain[i*ExcNum+k]:
                        Time = T/ms - MeanTime
                        if(abs(Time)<10):  LayerTemp[i][k].append(Time)
    PatternData.append(LayerTemp)
 
#RandomData = deepcopy(PatternData)
#for ii in range(PatternNum):
#    for i in range(LayerNum):
#        N = len(PatternData[ii][i])
#        for j in range(N):
#            Randii = random.randint(PatternNum)
#            RN = len(PatternData[Randii][i])
#            Randj = random.randint(RN)
#            Temp = RandomData[ii][i][j]
#            RandomData[ii][i][j] = RandomData[Randii][i][Randj]
#            RandomData[Randii][i][Randj] = Temp

 
if(RepTimes > 1):   
    MI = []
    LayerData = []
    for ii in range(LayerNum): 
        LayerDataNum = TestNum*LayerNum
        ThisLayer = []
        ILayer = ii
        for i in range(PatternNum):
            for j in range(RepTimes):
                ThisLayer.extend(ExcData[(i*RepTimes+j)*LayerDataNum +ILayer*TestNum \
                :(i*RepTimes+j)*LayerDataNum + (ILayer+1)*TestNum])                
        LayerData.append(ThisLayer)   
        TM = 0     
        q = 0.25
        for qq in range(1):
            CM = mi.TestInformation(ThisLayer,PatternNum,RepTimes,TestNum,q)
            if(CM > TM): TM = CM
            q *= 0.1
        MI.append(TM)   
    for ii in range(LayerNum):
        print("Information in layer %d is: %f"%((ii+1),MI[ii]))   
        
 
 
 
MaxPattern = 0            
MaxRep = 0
TheSignal = []
TheNoise = []

Signal = []
Noise = []
for i in range(LayerNum):
    Signal.append([])
    Noise.append([])               
for k in range(ExcNum):
    
    LayerSignal = []
    LayerNoise = []
    LayerAll = []
    for i in range(LayerNum):
        LayerSignal.append([])
        LayerNoise.append([])  
        LayerAll.append([])
        
    for i in range(LayerNum):
        for j in range(PatternNum):
            if(len(PatternData[j][i][k])>0):
                LayerNoise[i].append(std(PatternData[j][i][k]))
                LayerSignal[i].append(mean(PatternData[j][i][k]))
                
    tempPattern = 10000000
    tempRep = 10000000          
    for i in range(LayerNum):
        if(len(LayerSignal[i])>0):
            if (tempPattern > len(LayerSignal[i])): 
                tempPattern = len(LayerSignal[i])
            Signal[i].append(std(LayerSignal[i]))
        if(len(LayerNoise[i])>0):   
            if(tempRep > len(PatternData[0][i][k])): 
                tempRep = len(PatternData[0][i][k])
            Noise[i].append(mean(LayerNoise[i]))
            
    if(tempPattern > MaxPattern):
            MaxPattern = tempPattern
            TheSignal = copy.deepcopy(LayerSignal)
    if(tempRep > MaxRep):
            MaxRep = tempRep
            TheNoise = []
            for i in range(LayerNum):
                TheNoise.append(copy.deepcopy(PatternData[0][i][k]))           

for i in range(LayerNum):
    print("Layer%d: Signal=%f Noise=%f"%((i+1),mean(Signal[i]),mean(Noise[i])))   
    

#now = datetime.datetime.now()
#Datetime = now.strftime('%Y-%m-%d(%H-%M)') 
#
#if PatternNum > 1:
#    np.savetxt(Datetime+"FeedSignal.txt",mean(Signal,1))
#    for i in range(LayerNum):
#        np.savetxt(Datetime+"FeedPatternSpikesLayer%d.txt"%(i+1),TheSignal[i]) 
#if RepTimes > 1:
#    np.savetxt(Datetime+"FeedNoise.txt",mean(Noise,1))
#    for i in range(LayerNum):
#        np.savetxt(Datetime+"FeedNoiseSpikesLayer%d.txt"%(i+1),TheNoise[i]) 
#      

  
#%% Plots
#Neuron Position

subplot(3,1,1)
for i in range(LayerNum):
    plot(ExcU.t/ms, ExcU.v[i*ExcNum]/mV)
    plot(ExcU.t/ms, ExcU.v[i*ExcNum + 1]/mV,)  
xlim([0,SimulationTime/ms])      
subplot(3,1,2)
for i in range(LayerNum):
    plot(ExcU.t/ms, ExcU.gExc[i*ExcNum]/nS,)
    plot(ExcU.t/ms, ExcU.gInh[i*ExcNum]/nS,'--')
xlim([0,SimulationTime/ms])    
subplot(3,1,3)
plot(SucSpk.t/ms,SucSpk.i,'.b')    
plot(ExcSpk.t/ms,ExcSpk.i+200,'.k')   
plot(LatSpk.t/ms,LatSpk.i,'.r')   
xlim([0,SimulationTime/ms])

show()
