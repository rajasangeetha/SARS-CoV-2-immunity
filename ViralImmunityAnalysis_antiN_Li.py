##### 1. Import and get header data

import pandas as pd
import numpy

inputFile = "updated_OD_filev7_org.csv"

df2 = pd.read_csv(inputFile)

print(df2)

headers = list(df2.head(0))
headers

##### Take dataset without headers

noheaderdataset0 = pd.read_csv(inputFile).iloc[0:]
noheaderdataset0

##### Remove Rowws only with NaN => "NA:"

noheaderdataset = noheaderdataset0.dropna(axis = 0, how ='all')
noheaderdataset

print(len(noheaderdataset))

##### Take transpose of dataset

fulldataset = numpy.transpose(noheaderdataset) 

fulldataset

##### Get Length 
print(len(fulldataset))


####SARS-CoV-1 #####SARS-CoV-1 Average peak-normalized OD post-infection (3 mo.) 
#####Input Data: Average peak-normalized ELIZA ODs for SARS-CoV-1 N IgG antibodies from Li et al. 2006, "
#####Long-Term Persistence of Robust Antibody and Cytotoxic T Cell Responses in Recovered Patients Infected with SARS Coronavirus", PLoS ONE.
#####days : 30, 90, 180, 360, 540 #####ODs : 1.713299344, 1.636578538, 1.329728934, 0.767986943, 0.493304977

sarscov1data = [[1, 30, "SARSCoV1", 1.713299344, True, "NA", "NA"], [2, 90, "SARSCoV1", 1.238, False, (90 - 30), ((1.636578538-1.713299344) / 1.713299344)/(90-30)], [3, 180, "SARSCoV1", 1.16, False, (180 - 90), ((1.329728934-1.636578538) / 1.713299344)/(180 - 90)], [4, 360, "SARSCoV1", 1.043, False, (360 - 180), ((0.767986943-1.329728934) / 1.713299344)/(360 - 180)], [5, 540, "SARSCoV1", 0.919,False, (540 - 360), ((0.493304977-0.767986943) / 1.713299344)/(540 - 360)]]

print(sarscov1data)

##### SARS-CoV-1 Waning of Antibody OD

sarscov1datalength = len(sarscov1data)
print(sarscov1datalength)

##### Table/ArrayList creation

for i in range(0,len(sarscov1data)):
    for j in range(0,len(sarscov1data[i])):
        print(sarscov1data[i][j])
        
print(sarscov1data[2][6])

aod = 0.925
sarscov1paddedmeanwaning = dict();
sarscov1datalength3 = 5
index = 0

while aod <= 1.7:       
    index  = 2
    valueList = list()
    
    while index <= sarscov1datalength3:    
        print("aod :" + str(aod))
        if ( (sarscov1data[index-2][3] >= aod and aod >= sarscov1data[index-1][3]) or 
             (sarscov1data[index-2][3] <= aod and aod <= sarscov1data[index-1][3]) ) :
            valueList.append(sarscov1data[index-1][6])
        #else:
         #   valueList.append("None")
        sarscov1paddedmeanwaning["{0:.3f}".format(aod)] = valueList
            
        index = index + 1    
    aod = aod + 0.05

print(len(sarscov1paddedmeanwaning))

for key in sarscov1paddedmeanwaning:
    print(str(key) + " " + str(sarscov1paddedmeanwaning[key]))
    
import numpy as np

aodList = list(sarscov1paddedmeanwaning.keys())
print(aodList)

arraod = np.array(aodList)
arraod = arraod.astype(float)

covdataList = list(sarscov1paddedmeanwaning.values())
print(covdataList)

arrcov = np.array(covdataList).squeeze()

from scipy.interpolate import interp1d
print(len(aodList))
print(len(covdataList))

y_interpolation = interp1d(arraod, arrcov)

interpolate_x_new = 1.0


print("Value of Y at x = {} is".format(interpolate_x_new),
      y_interpolation(interpolate_x_new))
      
      
sarscov1antibodytimecourse = list()

sarscov1antibodytimecourse.append(1.0)

print(str(sarscov1antibodytimecourse))      

print(sarscov1antibodytimecourse[0])

##### Populating a dictionary for the interpolate values

day = 0
print("Populating a dictionary for the interpolate values")

while sarscov1antibodytimecourse[day] >= 0.925:
    print("For day : " +  str(day))
    print("sarscov1antibodytimecourse[day] : " + str(sarscov1antibodytimecourse[day]))
    print("sarscov1paddedmeanwaning : " + str(sarscov1paddedmeanwaning))
    print("inerpolate of sarscov1paddedmeanwaning")
    print(y_interpolation(sarscov1antibodytimecourse[day]))
    print("Added value")
    sarscov1antibodytimecourse.append(sarscov1antibodytimecourse[day] + y_interpolation(sarscov1antibodytimecourse[day]))    
    day = day + 1
    
print(sarscov1antibodytimecourse)

from matplotlib import pyplot as plt
plt.plot(sarscov1antibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('SARS-CoV-1')
plt.show()

import matplotlib.pyplot as plt
from datetime import datetime

f = plt.figure()
plt.plot(sarscov1antibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('SARS-CoV-1')
plt.show()

f.savefig("SARS-CoV-1_AntibodyTimecourse" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

###### Export to Excel

df3 = pd.DataFrame(sarscov1antibodytimecourse, columns=["Antibody Time Course"])
df3.to_csv('SARS-CoV-1-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + '.csv', index=False)

sarscov1baseline = 0.1298892

print(sarscov1baseline)

import math

def calculatelambda(sarscov1baseline, lambdaValue, index):
    exponentValue = -lambdaValue * index
    result = sarscov1baseline + (1 - sarscov1baseline) * math.exp(exponentValue)
    return result

sarscov1lsfuncwfixedbaseline = 0
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(sarscov1antibodytimecourse))):
    print(sarscov1antibodytimecourse[index] )
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.002
    
    ##Added 11/28
    sarscov1lsfuncwfixedbaseline = (sarscov1antibodytimecourse[index]) - math.pow((calculatelambda(sarscov1baseline, lambdaValue, index)), 2)     
    print(sarscov1lsfuncwfixedbaseline)
    
    if sarscov1lsfuncwfixedbaseline < currentMinValue:
        currentMinValue = sarscov1lsfuncwfixedbaseline      

print(currentMinValue)

sarscov1halflife = np.log(2) / 0.001039397687226661 

print(sarscov1halflife)

##### Added as part of python code

plt.plot(sarscov1antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('SARS-CoV-1 with fixed baseline')
plt.show()

##### Added as part of python code

f2 = plt.figure()
plt.plot(sarscov1antibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('SARS-CoV-1')
plt.show()

f2.savefig("SARS-CoV-1_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lambdaForPlot = 0.0010394
result = 0
plotBaseLines = list()

for days in range(0, len(sarscov1antibodytimecourse)):
    exponentValue = -lambdaForPlot * days
    result = sarscov1baseline + (1 - sarscov1baseline) * math.exp(exponentValue)
    plotBaseLines.append(result)
    
print(plotBaseLines)

from matplotlib import pyplot as plt

plt.plot(sarscov1antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('SARS-CoV-1 with fixed baseline')
#plt.show()



plt.title('SARS-CoV-1')
plt.show()

##### Added as part of python code

f2 = plt.figure()
plt.plot(plotBaseLines)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('SARS-CoV-1')
plt.show()

f2.savefig("SARS-CoV-1_withFixedBaseline_2" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

print(len(sarscov1antibodytimecourse))

sarscov1antibodytimecourseplusexp = sarscov1antibodytimecourse
day = len(sarscov1antibodytimecourse)

print(day)

#while day < 4393:
#Small file 86 records
#while day < 88:
maxDays = 86
while day < maxDays:
    day = day + 1
    print("day : " + day)
    exponentValueBaseLine = -0.001039397687226661
    tempValue = sarscov1baseline + (sarscov1antibodytimecourseplusexp[day-1] - sarscov1baseline) * math.exp(exponentValueBaseLine)
    print('tempValue' + str(tempValue))
    sarscov1antibodytimecourseplusexp.append(tempValue)
    
print(sarscov1antibodytimecourseplusexp)

print(len(sarscov1antibodytimecourseplusexp))

import matplotlib.pyplot as plt2
plt2.plot(sarscov1antibodytimecourseplusexp)
#f = plt2.figure()
plt2.ylim(0,1.1)
plt2.xlabel("Day")
plt2.ylabel("Peak-normalized antibody OD")
plt2.title('SARS-CoV-1')

import matplotlib.pyplot as plt3

f = plt3.figure()
plt3.plot(sarscov1antibodytimecourseplusexp)

plt3.ylim(0,1.1)
plt3.xlabel("Day")
plt3.ylabel("Peak-normalized antibody OD")
plt3.title('SARS-CoV-1')
plt3.show()

f.savefig("SARS-CoV-1_nOD-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

import math

def sarscov1probinfgivenaod(aod):
    exponentValue = 5.248568 + (9.749887 * aod)
    result = 1 / (1 + math.exp(exponentValue))
    return result
    
sarscov1probinfgivenaodList = list()

currentValue = 0;

while( currentValue <= 1):
    currentValue = currentValue + 0.00625
    sarscov1probinfgivenaodList.append(sarscov1probinfgivenaod(currentValue))

print(sarscov1probinfgivenaodList)

plt.plot(sarscov1probinfgivenaodList)
plt.ylim(0,0.006)
plt.xlim(0.0, 1.0)
plt.xlabel("Peak - normalized IgG antibody OD")
plt.ylabel("Daily probability of infection")
plt.title('SARS-CoV-1')
plt.show()

import matplotlib.pyplot as plt4

f2 = plt4.figure()
plt4.plot(sarscov1probinfgivenaodList)

plt4.ylim(0,0.006)
plt4.xlim(0.0, 1.0)
plt4.xlabel("Peak - normalized IgG antibody OD")
plt4.ylabel("Daily probability of infection")
plt4.title('SARS-CoV-1')
plt4.show()

f2.savefig("SARS-CoV-1_PrInf-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-1 Probability of Infection | Antibody OD

import math

def sarscov1probinfgivenall(a, b, aod):
    exponentValue = (-a) + (-b) * aod
    result = 1 / (1 + math.exp(exponentValue))
    return result
    
def populatesarscov1probinfList(a, b):
    sarscov1probinfList = list()
    
    x = 0;

    while( x <= 1):
        x = x + 0.05
        sarscov1probinfList.append(sarscov1probinfgivenall(a, b, x))
    
    return sarscov1probinfList  

from matplotlib import pyplot as plt

sarscov1probinfList1 = populatesarscov1probinfList(-5.248568, -9.749887)

plt.plot(sarscov1probinfList1, color='red')
plt.ylim(0,0.007)
plt.xlim(0.0, 1.0)
plt.xlabel("Peak-normalized IgG antibody OD")

sarscov1probinfList2 = populatesarscov1probinfList(-5.252705, -9.424764)
sarscov1probinfList3 = populatesarscov1probinfList(-4.593771, -12.8393)
sarscov1probinfList4 = populatesarscov1probinfList(-4.6620521, -12.79206)
sarscov1probinfList5 = populatesarscov1probinfList(-5.243759, -9.540128)

plt.plot(sarscov1probinfList2, color = 'green')
plt.plot(sarscov1probinfList3, color = 'lightblue')
plt.plot(sarscov1probinfList4, color = 'purple')
plt.plot(sarscov1probinfList5, color = 'blue')

plt.ylabel("Daily probability of infection")
plt.title('SARS-CoV-1')
#plt.show()

plt.title('SARS-CoV-1')
plt.show()


#####**************
#####The alternate values for a and b above come from our results using different approaches toward building the molecular evolutionary tree of the coronaviruses and toward building the time tree of the coronaviruses (see Supplement).

#####This panel goes into a supplementary figure.
#####*************

import matplotlib.pyplot as plt5

f5 = plt5.figure()
plt5.plot(sarscov1probinfList1, color='red')
plt5.plot(sarscov1probinfList2, color = 'green')
plt5.plot(sarscov1probinfList3, color = 'lightblue')
plt5.plot(sarscov1probinfList4, color = 'purple')
plt5.plot(sarscov1probinfList5, color = 'blue')

plt5.ylim(0,0.007)
plt5.xlim(0.0, 1.0)

plt5.xlabel("Peak - normalized IgG antibody OD")
plt5.ylabel("Daily probability of infection")
plt5.title('SARS-CoV-1')
plt5.show()

f5.savefig("SARS-CoV-1_PrInfs-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-1 Probability of No Reinfection Time Course

sarscov1probnoreinfectiontimecourse = list()
day = 0
sarscov1probnoreinfectiontimecourse.append(1.0)

sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1-sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day+1])))

print(sarscov1probnoreinfectiontimecourse)

print(len(sarscov1antibodytimecourse))

#### Need to change days to 4393
while day < 85:
    
    if(day < len(sarscov1antibodytimecourse)):
        sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1 - sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day + 1])))
    else:
        sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1 - sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[1392])))
            
    day = day + 1   

print(sarscov1probnoreinfectiontimecourse)

plt.plot(sarscov1probnoreinfectiontimecourse, color='blue')
#Limit not working
#plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Probability of no reinfection")
plt.title('SARS-CoV-1')
plt.show()

import matplotlib.pyplot as plt6

f6 = plt6.figure()
plt6.plot(sarscov1probnoreinfectiontimecourse, color = 'blue')

#plt6.ylim(0,0.007)
#plt6.xlim(0.0, 1.0)

plt6.xlabel("Day")
plt6.ylabel("Probability of no reinfection")
plt6.title('SARS-CoV-1')
plt6.show()

f6.savefig("SARS-CoV-1_PnorInfTimecourse-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

sarscov1probreinfection = list()

day = 0

#while day < 4393:
while day < 85:
    day = day + 1  
        
    #if day < 1392:
    if day < 42:
        tempValue = sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day]) * sarscov1probnoreinfectiontimecourse[day]
    else:
        tempValue = sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[3])* sarscov1probnoreinfectiontimecourse[day]
    
    sarscov1probreinfection.append(tempValue)      

import matplotlib.pyplot as plt7
plt7.plot(sarscov1probreinfection)
plt7.ylim(0, 0.0002)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('SARS-CoV-1') 

import matplotlib.pyplot as plt7

f = plt7.figure()
plt7.plot(sarscov1antibodytimecourseplusexp)

plt7.plot(sarscov1probreinfection)
plt7.ylim(0, 0.0002)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('SARS-CoV-1') 
f.savefig("SARS-CoV-1_PrInf-by-Day_reinfection" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

#####SARS-CoV-2

##SARS-CoV-2

##Input Data: Average peak-normalized ELIZA ODs for SARS-CoV-2 N IgG antibodies from Gudbjartsson et al. 2020, "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine.

##days : 35, 48, 70, 94, 109
##OD : 0.411, 0.399, 0.397, 0.350, 0.355

sarscov2data = [[1, 35, "SARSCoV2", 0.411/0.411, True, "NA", "NA"],
                [2, 48, "SARSCoV2", 0.399/0.411, False, (48 - 35), ((0.399 - 0.411) / 0.411)/(48 - 35)],
                [3, 70, "SARSCoV2", 0.397/0.411, False, (70 - 48), ((0.397 - 0.399) / 0.411)/(70 - 48)], 
                [4, 94, "SARSCoV2", 0.350/0.411, False, (95 - 70), ((0.350 - 0.397) / 0.411)/(94 - 70)]]             


## SARS-CoV-2 Waning of Antibody OD

sarscov2datalength = len(sarscov2data)

for i in range(0,len(sarscov2data)):
    for j in range(0,len(sarscov2data[i])):
        print(sarscov2data[i][j])
        
aod = 0.84
sarscov2paddedmeanwaning = dict();
sarscov2datalength = 4
index = 0

while aod <= 1:        
    index  = 2
    valueList = list()    
    
    while index <= sarscov2datalength:    
        print("aod :" + str(aod))
        if ( (sarscov2data[index-2][3] >= aod and aod >= sarscov2data[index-1][3]) or 
             (sarscov2data[index-2][3] <= aod and aod <= sarscov2data[index-1][3]) ) :
            valueList.append(sarscov1data[index-1][6])
        #else:
         #   valueList.append("None")
        if (len(valueList) != 0 ):
            sarscov2paddedmeanwaning["{0:.3f}".format(aod)] = valueList            
        index = index + 1    
    aod = aod + 0.02

print(len(sarscov2paddedmeanwaning))

for key in sarscov2paddedmeanwaning:
    print(str(key) + " " + str(sarscov2paddedmeanwaning[key]))
   
import numpy as np

aod2List = list(sarscov2paddedmeanwaning.keys())
aod2List

arraod2 = np.array(aod2List)
arraod2 = arraod2.astype(float)

covdata2List = list(sarscov2paddedmeanwaning.values())
covdata2List

arrcov2 = np.array(covdata2List).squeeze()

from scipy.interpolate import interp1d
print(len(aod2List))
print(len(covdata2List))

y_interpolation = interp1d(arraod2, arrcov2)

sarscov2antibodytimecourse = list()

sarscov2antibodytimecourse.append(1.0)

day = 0

while sarscov2antibodytimecourse[day] >= 0.86:
    print(day)
    print(sarscov2antibodytimecourse[day])
    print(sarscov2antibodytimecourse[day])
    print("sarscov2paddedmeanwaning")
    print(sarscov2paddedmeanwaning)
    print("sarscov2antibodytimecourse[day]")
    print(sarscov2antibodytimecourse[day])
    print("inerpolate of sarscov2paddedmeanwaning")
    print(y_interpolation(sarscov2antibodytimecourse[day]))
    print("Added value")
    test = sarscov2antibodytimecourse[day] + y_interpolation(sarscov2antibodytimecourse[day])
    sarscov2antibodytimecourse.append(sarscov2antibodytimecourse[day] + y_interpolation(sarscov2antibodytimecourse[day]))
    print("sarscov2antibodytimecourse")
    print(sarscov2antibodytimecourse)
    day = day + 1
    
    
print(sarscov2antibodytimecourse)

from matplotlib import pyplot as plt8
plt8.plot(sarscov2antibodytimecourse)
plt8.ylim(0,1.1)
plt8.xlabel("Day")
plt8.ylabel("Peak-normalized antibody OD")
plt8.title('SARS-CoV-1')
plt8.show()

df4 = pd.DataFrame(sarscov2antibodytimecourse, columns=["Antibody Time Course"])
df4.to_csv('SARS-CoV-2-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + '.csv', index=False)

sarscov2baseline = 0.1298892

###This baseline peak-normalized N IgG antibody level for SARS-CoV-2 comes from the ancestral 
#and descendent states analysis that used the baselines for the human "seasonal" coronaviruses 
#to estimate the baselines for the zoonotic coronaviruses.

sarscov2lsfuncwfixedbaseline = 0
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(sarscov2antibodytimecourse))):
    print(sarscov2antibodytimecourse[index] )
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.002
    ##Added 11/28
    sarscov2lsfuncwfixedbaseline = (sarscov2antibodytimecourse[index]) - math.pow((calculatelambda(sarscov2baseline, lambdaValue, index)), 2) 
    print(sarscov2lsfuncwfixedbaseline)
    
    if sarscov2lsfuncwfixedbaseline < currentMinValue:
        currentMinValue = sarscov2lsfuncwfixedbaseline      

print(currentMinValue)

sarscov2halflife = np.log(2) / 0.004661402990405273 

print(sarscov2halflife)

from matplotlib import pyplot as plt9
plt9.plot(sarscov2antibodytimecourse, color='red')
plt9.ylim(0,1.1)
plt9.xlabel("Day")
plt9.ylabel("Peak-normalized IgG antibody OD")
plt9.title('SARS-CoV-2')
plt9.show()

print(len(sarscov2antibodytimecourse))

sarscov2antibodytimecourseplusexp = sarscov2antibodytimecourse
day = len(sarscov2antibodytimecourse)

print(day)

#while day < 4393:
#Small file 86 records
#while day < 88:
maxDays = 86
while day < maxDays:
    day = day + 1
    print("day")
    print(day)
    exponentValueBaseLine = -0.004661402990405273
    tempValue = sarscov1baseline + (sarscov2antibodytimecourseplusexp[day-1] - sarscov1baseline) * math.exp(exponentValueBaseLine)
    print('tempValue')
    print(tempValue)
    sarscov2antibodytimecourseplusexp.append(tempValue)
    
print(sarscov2antibodytimecourseplusexp)

len(sarscov2antibodytimecourseplusexp)

import matplotlib.pyplot as plt2
plt2.plot(sarscov2antibodytimecourseplusexp)
#f = plt2.figure()
plt2.ylim(0,1.1)
plt2.xlabel("Day")
plt2.ylabel("Peak-normalized antibody OD")
plt2.title('SARS-CoV-2')

import matplotlib.pyplot as plt10

f = plt10.figure()
plt10.plot(sarscov2antibodytimecourseplusexp)

plt10.ylim(0,1.1)
plt10.xlabel("Day")
plt10.ylabel("Peak-normalized antibody OD")
plt10.title('SARS-CoV-1')
plt10.show()

f.savefig("SARS-CoV-2_nOD-by-Day_Antibodytimecours" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')  

##### SARS-CoV-2 Probability of Infection | Antibody OD

import math

def sarscov2probinfgivenaod(aod):
    exponentValue = 4.997881 + (11.099285 * aod)
    result = 1 / (1 + math.exp(exponentValue))
    return result
    
sarscov2probinfgivenaodList = list()

currentValue = 0;

while( currentValue <= 1):
    currentValue = currentValue + 0.00625
    sarscov2probinfgivenaodList.append(sarscov2probinfgivenaod(currentValue))
    
print(sarscov2probinfgivenaodList)

import matplotlib.pyplot as plt11

plt11.plot(sarscov2probinfgivenaodList)
plt11.ylim(0,0.007)
plt11.xlim(0.0, 1.0)
plt11.xlabel("Peak - normalized IgG antibody OD")
plt11.ylabel("Daily probability of infection")
plt11.title('SARS-CoV-2')
plt11.show()

import matplotlib.pyplot as plt12

f12 = plt12.figure()
plt12.plot(sarscov2probinfgivenaodList)

plt12.ylim(0,0.007)
plt12.xlim(0.0, 1.0)
plt12.xlabel("Peak - normalized IgG antibody OD")
plt12.ylabel("Daily probability of infection")
plt12.title('SARS-CoV-2')
plt12.show()

f12.savefig("SARS-CoV-2_PrInf-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

from matplotlib import pyplot as plt13

sarscov2probinfList1 = populatesarscov1probinfList(-4.997881, -11.099285)

plt13.plot(sarscov1probinfList1, color='red')
plt13.ylim(0,0.008)
plt13.xlim(0.0, 1.0)
plt13.xlabel("Peak-normalized IgG antibody OD")

sarscov2probinfList2 = populatesarscov1probinfList(-5.1806, -9.802874)
sarscov2probinfList3 = populatesarscov1probinfList(-4.76283, -12.00514)
sarscov2probinfList4 = populatesarscov1probinfList(-4.5831338, -13.20581)
sarscov2probinfList5 = populatesarscov1probinfList(-5.12873, -10.142186)

plt13.plot(sarscov2probinfList2, color = 'green')
plt13.plot(sarscov2probinfList3, color = 'lightblue')
plt13.plot(sarscov2probinfList4, color = 'purple')
plt13.plot(sarscov2probinfList5, color = 'blue')

plt13.ylabel("Daily probability of infection")
plt13.title('SARS-CoV-2')
#plt.show()

import matplotlib.pyplot as plt14

f14 = plt14.figure()
plt14.plot(sarscov2probinfList1, color='red')
plt14.plot(sarscov2probinfList2, color = 'green')
plt14.plot(sarscov2probinfList3, color = 'lightblue')
plt14.plot(sarscov2probinfList4, color = 'purple')
plt14.plot(sarscov2probinfList5, color = 'blue')

plt14.ylim(0,0.007)
plt14.xlim(0.0, 1.0)

plt14.xlabel("Peak - normalized IgG antibody OD")
plt14.ylabel("Daily probability of infection")
plt14.title('SARS-CoV-1')
plt14.show()

f14.savefig("SARS-CoV-2_PrInfs-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-2 Probability of No Reinfection Time Course

sarscov2probnoreinfectiontimecourse = list()

day = 0
sarscov2probnoreinfectiontimecourse.append(1.0)

sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1-sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day+1])))

print(sarscov2probnoreinfectiontimecourse)

print(len(sarscov2antibodytimecourse))

#### Need to change days to 4393
while day < 85:
    
    if(day < len(sarscov2antibodytimecourse)):
        sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day + 1])))
    else:
        #sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1392])))
        sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[42])))
            
    day = day + 1   
    
print(sarscov2probnoreinfectiontimecourse)

plt.plot(sarscov2probnoreinfectiontimecourse, color='blue')
#Limit not working
#plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Probability of no reinfection")
plt.title('SARS-CoV-2')
plt.show()

import matplotlib.pyplot as plt15

f15 = plt15.figure()
plt15.plot(sarscov2probnoreinfectiontimecourse, color = 'blue')

#plt6.ylim(0,0.007)
#plt6.xlim(0.0, 1.0)

plt15.xlabel("Day")
plt15.ylabel("Probability of no reinfection")
plt15.title('SARS-CoV-1')
plt15.show()

f15.savefig("SARS-CoV-2_PnorInfTimecourse-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

sarscov2probreinfection = list()
day = 0;

#while day < 4393:
while day < 85:
    day = day + 1  
        
    #if day < 1392:
    if day < 42:
        tempValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day]) * sarscov2probnoreinfectiontimecourse[day]
    else:       
        #empValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1392])  * sarscov2probnoreinfectiontimecourse[day]
        tempValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[43])  * sarscov2probnoreinfectiontimecourse[day]
    
    sarscov2probreinfection.append(tempValue)      
    
print(sarscov2probreinfection)

import matplotlib.pyplot as plt17
plt17.plot(sarscov2probreinfection)
plt17.ylim(0, 0.0005)
plt17.xlabel("Day")
plt17.ylabel("Probability of reinfection")
plt17.title('SARS-CoV-2') 

import matplotlib.pyplot as plt17

f17 = plt17.figure()
plt17.plot(sarscov2probreinfection)

plt17.ylim(0,0.005)
plt17.xlim(0.0, 1.0)
plt17.xlabel("Peak - normalized IgG antibody OD")
plt17.ylabel("Daily probability of infection")
plt17.title('SARS-CoV-2')
plt17.show()

f17.savefig("SARS-CoV-2_PrInf-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')