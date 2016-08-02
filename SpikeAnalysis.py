# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 20:51:36 2015

@author: Hong
"""
import numpy as np
from copy import deepcopy




#%%------------------------Information Analysis --------------------------%%#
def Distance(Sa,Sb,q):
    m = len(Sa)
    n = len(Sb)
    Sheet = np.zeros([m+1,n+1])
    Sheet[:,0] = range(m+1)
    Sheet[0,:] = range(n+1)
    
    for i in range(m):
        for j in range(n):
            G1 = Sheet[i,j+1] + 1
            G2 = Sheet[i+1,j] + 1
            G3 = Sheet[i,j] + q*abs(Sa[i] - Sb[j])
            Sheet[i+1,j+1] = min([G1,G2,G3])
    Value = Sheet[m,n]        
    return Value


def GroupDistance(GS, TS, NeuronNum, q): #GroupSquence, TargetSquence  
    MemberNum = len(GS)
    Ds = []
    for i in range(NeuronNum):
        Dist = 0
        for j in range(MemberNum):           
            if(len(TS[i])+len(GS[j][i])>0):
                Dist += (Distance(GS[j][i],TS[i],q)+0.1)**-0.1
#            except IndexError:
#                print("List Index i=%d, j=%d out of range"%(i,j))      
        Ds.append((Dist)**10)#-1/2

            
    MeanDistance = np.mean(Ds)    
    return MeanDistance**-2#-2
    
    
def Clustering(Groups, TS, NeuronNum, q, softDist=0):
    ''' Dscription
    '''
    GroupNum = len(Groups)
    Ds = []
    for i in range(GroupNum):
        Ds.append(GroupDistance(Groups[i],TS,NeuronNum,q))
    MinDist = min(Ds)
    MeanDist = np.mean(Ds)
    DistPower = MeanDist-MinDist
    Ds = np.array(Ds)
    Nc = np.zeros(GroupNum)
    if(softDist==0):
        Ind = np.equal(Ds,MinDist)
        k = np.count_nonzero(Ind)
        Nc[Ind] = 1/float(k)
    elif softDist > 0 :
        Ind = Ds<=MinDist*(1+softDist)
        k = sum(1/(Ds[Ind]+1e-5))
        Nc[Ind] = 1/((Ds[Ind]+1e-5)*k)
    return [Nc,DistPower]    
        
def ConfusingMatrix(SpikeTrains, StimulusNum, RepeatNum, NeuronNum, q, softDist=0):
    ''' Compute Clustering Distance ...
        ''' 
    TrailNum = len(SpikeTrains)
    if TrailNum != StimulusNum*RepeatNum*NeuronNum:
        return False
    #---------------Grouping----------------#
    Groups = []
    Trail = 0
    for i in range(StimulusNum):
        GroupTrains = []
        for j in range(RepeatNum):
            Population = []
            for k in range(NeuronNum):
                Population.append(SpikeTrains[Trail])
                Trail+=1
            GroupTrains.append(Population)
        Groups.append(GroupTrains)
#    Quite = []
#    for i in range(NeuronNum):
#         Quite.append([]) 
#    Groups.append([Quite])
   #---------------Clustering---------------#
    Matrix =  np.zeros([StimulusNum,StimulusNum])
    Dist = []
    for i in range(StimulusNum):
        for j in range(RepeatNum):
            #print("runing the %dth Stimulus and %dth times"%(i,j))
            TempGroup = deepcopy(Groups)      
            for kk in range(StimulusNum):
                del TempGroup[kk][j]
            [Nc,DistPower] = Clustering(TempGroup,Groups[i][j],NeuronNum,q,softDist)
            Dist.append(DistPower)
            Matrix[i,:] += Nc
    
#    Matrix[StimulusNum,StimulusNum] = RepeatNum         
   #--------------Return-------------------#
   # print("Dist Different=%f"%(np.mean(Dist)))
    return Matrix
   
def log(x):
    if x<=0:
        return 1e10
    else:
        return np.log2(x)
            
        
   
   
def TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum, q, softDist=0):
    ''' ... '''
    #Native Information
    Mat = ConfusingMatrix(SpikeTrains,StimulusNum,RepeatNum,NeuronNum,q,softDist)
    H = 0
    totNum = (StimulusNum)*RepeatNum
    for i in range(StimulusNum):
        for j in range(StimulusNum):
            H += Mat[i,j]*(log(Mat[i,j])-log(sum(Mat[:,j])) \
                     -log(sum(Mat[i,:]))+log(totNum))
       
    H = H/totNum
    
#    #Information Bias 
#    Length = len(SpikeTrains)
#    Indx = np.random.permutation(Length)
#    MassTrains = []
#    for i in range(Length):
#        MassTrains.append(SpikeTrains[Indx[i]])
#    Mat = ConfusingMatrix(MassTrains,StimulusNum,RepeatNum,NeuronNum,q,softDist)
#    Hb = 0
#    totNum = (StimulusNum)*RepeatNum
#    for i in range(StimulusNum):
#        for j in range(StimulusNum):
#            Hb += Mat[i,j]*(log(Mat[i,j])-log(sum(Mat[:,j])) \
#                     -log(sum(Mat[i,:]))+log(totNum))
#       
#    Hb = Hb/totNum
    
    #print("H = %f, Hb=%f"%(H,Hb))
    print("H = %f"%(H))
    return H
def BiasInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum, q, softDist=0):   
    #Information Bias 

    Length = len(SpikeTrains)
    Indx = np.random.permutation(Length)
    MassTrains = []
    for i in range(Length):
        MassTrains.append(SpikeTrains[Indx[i]])
    Mat = ConfusingMatrix(MassTrains,StimulusNum,RepeatNum,NeuronNum,q,softDist)
    Hb = 0
    totNum = (StimulusNum)*RepeatNum
    for i in range(StimulusNum):
        for j in range(StimulusNum):
            Hb += Mat[i,j]*(log(Mat[i,j])-log(sum(Mat[:,j])) \
                     -log(sum(Mat[i,:]))+log(totNum))
       
    Hb = Hb/totNum
    print("Hb = %f"%(Hb))    
    return Hb
    
              
def Information(SpikeTrains, StimulusNum, RepeatNum, NeuronNum, softDist=0):
    kk = 0.38196601125010515179541316563436
    OL = 0.001
    OR = 5
    QL = OL + (OR-OL)*kk
    QR = OR - (OR-OL)*kk
    
    OHL = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,OL, softDist)
    OHR = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,OR, softDist)
    HL = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QL, softDist)
    HR = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QR, softDist)
    for i in range(200):
        print("Test %d times QL=%f and QR=%f"%(i+1,QL,QR))
        if  abs(OR-OL) <0.5 or (HR==HL and HR>=OHR and HL>=OHL):
            break
        elif HR> OHR:
             OR = QR
             QR = QL
             QL = OL + kk*(OR-OL)
             OHR =HR
             HR = HL
             HL = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QL, softDist)
        elif HL>OHL:
             OL = QL
             QL = QR
             QR = OR - kk*(OR-OL)
             OHL = HL
             HL = HR
             HR = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QR, softDist)
        elif OHL>OHR:
             OR = QR
             QR = QL
             QL = OL + kk*(OR-OL)
             OHR =HR
             HL = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QL, softDist)
        else:
             OL = QL
             QL = QR
             QR = OR - kk*(OR-OL)
             OHL = HL
             HL = HR
             HR = TestInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QR, softDist)            
                
    MAXH = np.max([HR,HL,OHR,OHL])            
    if HR == MAXH:
        #HB = BiasInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QR, softDist)   
        return [HR, QR]
    elif HL == MAXH:
        #HB = BiasInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,QL, softDist) 
        return [HL, QL]
    elif OHL == MAXH:
        #HB = BiasInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,OL, softDist) 
        return [OHL, OL]
    else:
        #HB = BiasInformation(SpikeTrains, StimulusNum, RepeatNum, NeuronNum,OR, softDist) 
        return [OHR, OR]
   
   
#%%    
#ExcFile= np.load("ExcData_2016-01-06(11-36).npz") 
#Data = ExcFile['arr_0']
#
#NeuronNum=400
#LayerNum =4
#PatternNum=8
#RepTimes= 20
#
#LayerDataNum = NeuronNum*LayerNum
#LayerData = []
#ILayer = 1
#for i in range(PatternNum):
#    for j in range(RepTimes):
#        LayerData.extend(Data[(i*RepTimes+j)*LayerDataNum +ILayer*NeuronNum \
#        :(i*RepTimes+j)*LayerDataNum + (ILayer+1)*NeuronNum])
#        
##M = ConfusingMatrix(LayerData,PatternNum,RepTimes,NeuronNum,0.5,0)   
#H = Information(LayerData,PatternNum,RepTimes,NeuronNum,0)
##print("information=%f, q=%f"%(H[0],H[1]))           
#   
   
   
   
#%%--------------------------Spike Pattern Statsics-------------------------%%#  
def SpikeTrainSort(Spikes):
    Num = len(Spikes)
    
    SpikeTimes = []
    for i in range(Num):
        for T in Spikes[i]:
            SpikeTimes.append(T*1000)
  
    SpikeTimes = np.sort(SpikeTimes)
    meanTime = 0
    Num = len(SpikeTimes)
    halfNum = int(Num/2)
    meanTime=SpikeTimes[halfNum]
    FilteredTimes = []
    Num = len(SpikeTimes)
    for T in SpikeTimes:
        if abs(T - meanTime) < 5:
            FilteredTimes.append(T)
    return FilteredTimes   

def VolleyStatistic(Spikes):
    SpikeTimes = SpikeTrainSort(Spikes)
    sigma = np.std(SpikeTimes)
    num = len(SpikeTimes)
    return [num,sigma]
    
    
   
   
   
    
       
