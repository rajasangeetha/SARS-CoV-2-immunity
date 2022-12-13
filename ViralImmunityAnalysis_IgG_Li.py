##### Durability of Immunity against reinfection of SARS-CoV-2

##### Author: Sangeetha Vijayam
##### Reviwer: Jeffrey P. Townsend, Alex Dornburg, Hayley B. Hassler
##### Date: 12/01/2022

##### SARS-CoV-1


##### Input Data: updated_OD _filev7.csv contains processed N IgG antibodies converted to IgG antibodies from the supplementary dataset 
#of Edridge et al. 2020, "Seasonal coronavirus protective immunity is short-lasting", Nature Medicine. using a 
#linear regression model of N and IgG antibody data from Li et al. 2006, "Long-Term Persistence of Robust Antibody and 
#Cytotoxic T Cell Responses in Recovered Patients Infected with SARS Coronavirus", PLoS ONE.

##### 1. Import and get header data

import pandas as pd
import numpy
from datetime import datetime
import math

inputFile = "updated_OD_filev7_org.csv"

df2 = pd.read_csv(inputFile)

print(df2)

headers = list(df2.head(0))
print(headers)

##### Take dataset without headers

noheaderdataset0 = pd.read_csv(inputFile).iloc[0:]
noheaderdataset0

##### Remove Rowws only with NaN => "NA:"

noheaderdataset = noheaderdataset0.dropna(axis = 0, how ='all')
noheaderdataset

print(len(noheaderdataset))

##### Take transpose of dataset

print(noheaderdataset)

edridgefulldataset = numpy.transpose(noheaderdataset) 

edridgefulldataset 

edridgefulldataset = numpy.transpose(edridgefulldataset)

print(edridgefulldataset)

##### Get Length 

print(len(edridgefulldataset))

#####  229E
###### Input Data: ELIZA ODs for HCoV-229E N IgG antibodies converted to IgG antibodies from the supplementary 
###### dataset of Edridge et al. 2020, "Seasonal coronavirus protective immunity is short-lasting", Nature Medicine. 
###### using a linear regression model of N and IgG antibody data from Li et al. 2006, "Long-Term Persistence of 
###### Robust Antibody and Cytotoxic T Cell Responses in Recovered Patients Infected with SARS Coronavirus", PLoS ONE.

x229edata = df2[df2.isin(["229E"]).any(axis=1)]

print(x229edata.iloc[:1])

llx229elength = len(x229edata)

###### Note: Have to check Full Simplify

##### Probability of Infection | Antibody OD

###### These values of a and b go into the the ancestral and descendent states analysis to estimate a and b for the zoonotic coronaviruses.

#Yet to update
a = 1
b =2

import math

def x229eprobinfgivenaod(aod, a, b):
    exponentValue = -(a + b * aod)   
    result = 1 / (1 + math.exp(exponentValue))
    return result
	
x229probinfgivenaodList = list()

currentValue = 0;

while( currentValue <= 1):
    currentValue = currentValue + 0.00625
    x229probinfgivenaodList.append(x229eprobinfgivenaod(currentValue, a, b))

print(x229probinfgivenaodList)


import matplotlib.pyplot as plt

plt.plot(x229probinfgivenaodList)
plt.ylim(0,0.006)
plt.xlim(0.0, 1.0)
plt.xlabel("Peak - normalized IgG antibody OD")
plt.ylabel("Daily probability of infection")
plt.title('SARS-CoV-1')
plt.show()

import matplotlib.pyplot as plt4

f2 = plt4.figure()
plt4.plot(x229probinfgivenaodList)

plt4.ylim(0,0.006)
plt4.xlim(0.0, 1.0)
plt4.xlabel("Peak - normalized IgG antibody OD")
plt4.ylabel("Daily probability of infection")
plt4.title('SARS-CoV-1')
plt4.show()

f2.savefig("SARS-CoV-1_PrInf-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

print(x229edata)

x229edataArr = x229edata.to_numpy()

import statistics
import numpy as np

def grade_avg(num_list):
   
    grades = [float(num) for num in num_list]
   
    try:
        return statistics.mean(grades)
    except:
        return np.nan
		
import re

aod = 0
x229epaddedmeanwaning = dict();
llx229elength = len(x229edata)

index = 0

while aod <= 2.9:       
    index  = 2
    valueList = list()
    
    while index <= llx229elength:    
        #print("aod :" + str(aod))
        #print("index : " + str(index))
        print(x229edataArr[index-2])
        
        if ((x229edataArr[index-2][9] >= aod and aod >= x229edataArr[index-2][9]) or 
            (x229edataArr[index-2][9] <= aod and aod <= x229edataArr[index-2][9]) ) :
            if (x229edataArr[[index-1, 10]] != "NA" and re.match( r'^(F).*(E)$', x229edataArr[index-2][6] )):
                if x229edataArr[index-1][8].isnumeric() == True: 
                    valueList.append(x229edataArr[index-1][8])
                else:
                    valueList.append(0)
        average = grade_avg(valueList)
        print("average : " + str(average))
        x229epaddedmeanwaning["{0:.3f}".format(aod)] = average
        
        index = index + 1    
    aod = aod + 0.05

print(x229epaddedmeanwaning)

for i in range(0, 2):
    (k := next(iter(x229epaddedmeanwaning)), x229epaddedmeanwaning.pop(k))

print(x229epaddedmeanwaning)

for i in range(0, 2):
    (k := next(iter(x229epaddedmeanwaning)), x229epaddedmeanwaning.pop(k))

print(x229epaddedmeanwaning)

updict = {'0.05':'0', '0':'0'}

# ** operator for packing and unpacking items in order
x229epaddedmeanwaningNew = {**updict, **x229epaddedmeanwaning}
print(x229epaddedmeanwaningNew)

import numpy as np

aodList = list(x229epaddedmeanwaningNew.keys())
print(aodList)

arraod = np.array(aodList)
arraod = arraod.astype(float)

covdataList = list(x229epaddedmeanwaningNew.values())
print(covdataList)

arrcov = np.array(covdataList).squeeze()

from scipy.interpolate import interp1d
print(len(aodList))
print(len(covdataList))

y_interpolation = interp1d(arraod, arrcov)

interpolate_x_new = 1.0


print("Value of Y at x = {} is".format(interpolate_x_new),
      y_interpolation(interpolate_x_new))
      
x229eantibodytimecourse = list()

x229eantibodytimecourse.append(1.0)

print(str(x229eantibodytimecourse))

day = 0
print("Populating a dictionary for the interpolate values")

while x229eantibodytimecourse[day] >= 0:
    print("For day : " +  str(day))
    print("x229eantibodytimecourse[day] : " + str(x229eantibodytimecourse[day]))
    print("x229eantibodytimecourse : " + str(x229eantibodytimecourse))
    print("inerpolate of x229epaddedmeanwaningNew")
    print(y_interpolation(x229eantibodytimecourse[day]))
    print("Added value")
    x229eantibodytimecourse.append(x229eantibodytimecourse[day] + y_interpolation(x229eantibodytimecourse[day]))    
    day = day + 1
    
print(x229eantibodytimecourse)

from matplotlib import pyplot as plt
plt.plot(x229eantibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
#plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-229E')
plt.show()

import matplotlib.pyplot as plt
import math

f = plt.figure()
plt.plot(x229eantibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
#plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-229E')
plt.show()

f.savefig("HCoV-229E_AntibodyTimecourse" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

import math

def calculatelambda(sarscov1baseline, lambdaValue, index):
    exponentValue = -lambdaValue * index
    result = sarscov1baseline + (1 - sarscov1baseline) * math.exp(exponentValue)
    return result
    
x229elsfunc = 0
baseline = 0.1
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(x229eantibodytimecourse))):
    print(x229eantibodytimecourse[index])
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.005
    
    ##Added 11/28
    nl63lsfunc = (x229eantibodytimecourse[index]) - math.pow((calculatelambda(baseline, lambdaValue, index)), 2)     
    print(x229elsfunc)
    
    if x229elsfunc < currentMinValue:
        currentMinValue = x229elsfunc    
        
print(currentMinValue)

##### This estimate of lambda for HCoV-229E of 0.0999476 goes into the the ancestral and descendent states analysis.

x229ehalflife  = np.log(2) / 0.09994762543639751

print(x229ehalflife)

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.0)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E')
plt.show()

f2 = plt.figure()
plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E')
plt.show()

f2.savefig("HCoV-229E_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lambdaForPlot = 0.00231669
baseline = 0.249782
result = 0
plotBaseLines = list()

for days in range(0, 1675):
    exponentValue = -lambdaForPlot * days
    result = baseline + (1 - baseline) * math.exp(exponentValue)
    plotBaseLines.append(result)

print(plotBaseLines)

from matplotlib import pyplot as plt

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-NL63 with fixed baseline')
#plt.show()

plt.title('HCoV-229E')
plt.show()

f2 = plt.figure()
plt.plot(plotBaseLines)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-229E')
plt.show()

f2.savefig("HCoV-229E_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

from datetime import datetime

df3 = pd.DataFrame(x229eantibodytimecourse, columns=["Antibody Time Course"])
df3.to_csv('HCoV-229E-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + '.csv', index=False)

print(len(x229eantibodytimecourse))

x229elsfunc = 0
baseline = 0.0999476
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(x229eantibodytimecourse))):
    print(x229eantibodytimecourse[index])
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.00631483
    
    ##Added 11/28
    x229elsfunc = (x229eantibodytimecourse[index]) - math.pow((calculatelambda(baseline, lambdaValue, index)), 2)     
    print(x229elsfunc)
    
    if x229elsfunc < currentMinValue:
        currentMinValue = x229elsfunc      
        
print(currentMinValue)

##### This estimate of lambda for HCoV-229E of 0.0999476 goes into the the ancestral and descendent states analysis.

x229ehalflife = np.log(2) / 0.006314829758653257

print(x229ehalflife)

##### Added as part of python code

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.0)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E')
plt.show()

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.0)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E')
plt.show()

f2 = plt.figure()
plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E')
plt.show()

f2.savefig("HCoV-229E_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lambdaForPlot = 0.00631483
baseline = 0.0999476
result = 0
plotBaseLines = list()

for days in range(0, len(x229eantibodytimecourse)):
    exponentValue = -lambdaForPlot * days
    result = baseline + (1 - baseline) * math.exp(exponentValue)
    plotBaseLines.append(result)
    
print(plotBaseLines)

from matplotlib import pyplot as plt

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-OC43 with fixed baseline')
#plt.show()

plt.title('HCoV-OC43')
plt.show()


f2 = plt.figure()
plt.plot(plotBaseLines)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-OC43')
plt.show()

f2.savefig("HCoV-OC43_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

print(len(x229eantibodytimecourse))

x229eantibodytimecourseplusexp = x229eantibodytimecourse
day = len(x229eantibodytimecourse)

print(day)

#while day < 4393 - Length[oc43antibodytimecourse]:
#Small file 86 records
#while day < 89:
maxDays = 2
while day < maxDays:
    day = day + 1
    print("day : " + str(day))
    exponentValueBaseLine = -0.006314829758653257
    tempValue = 0.09995743181844773 + (x229eantibodytimecourseplusexp[day-1] - 0.09994762543639751) * math.exp(exponentValueBaseLine)
    print('tempValue' + str(tempValue))
    x229eantibodytimecourseplusexp.append(tempValue)
    
print(x229eantibodytimecourseplusexp)

import matplotlib.pyplot as plt2
plt2.plot(x229eantibodytimecourseplusexp)
plt2.ylim(0,1.1)
plt2.xlabel("Day")
plt2.ylabel("Peak-normalized antibody OD")
plt2.title('HCoV-229E')

import matplotlib.pyplot as plt3

f = plt3.figure()
plt3.plot(x229eantibodytimecourseplusexp)

plt3.ylim(0,1.1)
plt3.xlabel("Day")
plt3.ylabel("Peak-normalized antibody OD")
plt3.title('HCoV-OC43')
plt3.show()

f.savefig("HCoV-229E_nOD-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lengthabtc = len(x229eantibodytimecourse)
print(lengthabtc)

##### x229e Probability of No Reinfection Time Course

x229eprobnoreinfectiontimecourse = list()
day = 0
x229eprobnoreinfectiontimecourse.append(1.0)

x229eprobnoreinfectiontimecourse.append(x229eprobnoreinfectiontimecourse[day] * (1-x229eprobinfgivenaod(x229eantibodytimecourse[day], a, b)))

print(x229eprobnoreinfectiontimecourse)

#### Need to change days to 4393
while day < 1:
    
    if(day < lengthabtc):
        x229eprobnoreinfectiontimecourse.append(x229eprobnoreinfectiontimecourse[day] * (1 - x229eprobinfgivenaod(x229eantibodytimecourse[day + 1], a, b)))
    else:
        x229eprobnoreinfectiontimecourse.append(x229eprobnoreinfectiontimecourse[day] * (1 - x229eprobinfgivenaod(x229eantibodytimecourse[1392], a, b)))
            
    day = day + 1   
    
print(x229eprobnoreinfectiontimecourse)

plt.plot(x229eprobnoreinfectiontimecourse, color='blue')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Probability of no reinfection")
plt.title('HCoV-229E')
plt.show()

import matplotlib.pyplot as plt6

f6 = plt6.figure()
plt6.plot(x229eprobnoreinfectiontimecourse, color = 'blue')

#plt6.ylim(0,0.007)
#plt6.xlim(0.0, 1.0)

plt6.xlabel("Day")
plt6.ylabel("Probability of no reinfection")
plt6.title('HCoV-229E')
plt6.show()

f6.savefig("HCoV-229E_PnorInfTimecourse-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

x229eprobreinfection = list()

day = 0

#while day < 4393:
while day < 2:
    day = day + 1  
        
    #if day < 1392:
    if day < 2:
        tempValue = x229eprobinfgivenaod(x229eantibodytimecourse[day], a, b) * x229eantibodytimecourse[day]
    else:
        #tempValue = oc43probinfgivenaod(oc43antibodytimecourse[lengthabtc])* oc43probnoreinfectiontimecourse[day]
        tempValue = x229eprobinfgivenaod(x229eantibodytimecourse[1], a, b)* x229eantibodytimecourse[1]
    
    x229eprobreinfection.append(tempValue)      

import matplotlib.pyplot as plt7
plt7.plot(x229eprobreinfection)
plt7.ylim(0, 1.0)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('HCoV-229E') 

import matplotlib.pyplot as plt7

f = plt7.figure()

plt7.plot(x229eprobreinfection)
plt7.ylim(0, 1.0)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('HCoV-229E') 
f.savefig("HCoV-229E_PrInf-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-1

##### SARS-CoV-1 Average peak-normalized OD post-infection (3 mo.)

##### Input Data: Average peak-normalized ELIZA ODs for SARS-CoV-1 IgG antibody data from Li et al. 2006, 
#"Long-Term Persistence of Robust Antibody and Cytotoxic T Cell Responses in Recovered Patients Infected with 
#SARS Coronavirus", PLoS ONE.

##### days : 90, 180, 360, 540, 720 ODs : 0.95351472, 0.88161949, 0.71702398, 0.56000208, 0.51459265

sarscov1data = [[1, 90, "SARSCoV1", 0.95351472, True, "NA", "NA"], 
                [2, 180, "SARSCoV1", 0.88161949, False, (180 - 90), ((0.88161949 - 0.95351472) / 0.95351472) / (180 - 90)],
                [3, 360, "SARSCoV1", 0.71702398, False, (360 - 180), ((0.71702398- 0.88161949) / 0.95351472) / (360 - 180)],
                [4, 540, "SARSCoV1", 0.56000208, False, (540 - 360),  ((0.56000208 - 0.71702398) / 0.95351472) / (540 - 360)],
                [5, 720, "SARSCoV1", 0.51459265, False, (720 - 540),  ((0.51459265 - 0.56000208) / 0.95351472) / (720 - 540)]]
                
                
                
###### SARS-CoV-1 Waning of Antibody OD

sarscov1datalength = len(sarscov1data)

print(sarscov1datalength)

for i in range(0,len(sarscov1data)):
    for j in range(0,len(sarscov1data[i])):
        print(sarscov1data[i][j])
        
        

aod = 0.525
sarscov1paddedmeanwaning = dict();
sarscov1datalength = 4
index = 0

while aod <= 0.95:        
    index  = 2
    valueList = list()    
    
    while index <= sarscov1datalength:    
        print("aod :" + str(aod))
        if ( (sarscov1data[index-2][3] >= aod and aod >= sarscov1data[index-1][3]) or 
             (sarscov1data[index-2][3] <= aod and aod <= sarscov1data[index-1][3]) ) :
            valueList.append(sarscov1data[index-1][6])
        #else:
         #   valueList.append("None")
        if (len(valueList) != 0 ):
            sarscov1paddedmeanwaning["{0:.3f}".format(aod)] = valueList            
        index = index + 1    
    aod = aod + 0.025

print(len(sarscov1paddedmeanwaning))


sarscov1meanwaning = sarscov1paddedmeanwaning

for key in sarscov1paddedmeanwaning:
    print(str(key) + " " + str(sarscov1paddedmeanwaning[key]))
    
import numpy as np

aod2List = list(sarscov1paddedmeanwaning.keys())
aod2List

arraod2 = np.array(aod2List)
arraod2 = arraod2.astype(float)

covdata2List = list(sarscov1paddedmeanwaning.values())
covdata2List

arrcov2 = np.array(covdata2List).squeeze()

from scipy.interpolate import interp1d
print(len(aod2List))
print(len(covdata2List))

y_interpolation = interp1d(arraod2, arrcov2)

sarscov1antibodytimecourse = list()

sarscov1antibodytimecourse.append(0.95)

print(sarscov1antibodytimecourse)

day = 0

#while ( day <= (len(sarscov1antibodytimecourse) - 1)  and sarscov1antibodytimecourse[day] >= 0.525):
while ( sarscov1antibodytimecourse[day] <= 0.525):
    print(day)
    print("sarscov1antibodytimecourse[day] : " + str(sarscov1antibodytimecourse[day]))
    print("sarscov1paddedmeanwaning : " + str(sarscov1paddedmeanwaning))
    print("inerpolate of sarscov1paddedmeanwaning")
    #print(y_interpolation(sarscov1antibodytimecourse[day]))
    print("Added value")
    #sarscov1antibodytimecourse.append(sarscov1antibodytimecourse[day] + y_interpolation(sarscov1antibodytimecourse[day]))  
    sarscov1antibodytimecourse.append(sarscov1antibodytimecourse[day] + 0)
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
import math
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

from datetime import datetime
import pandas as pd

df3 = pd.DataFrame(sarscov1antibodytimecourse, columns=["Antibody Time Course"])
df3.to_csv('SARS-CoV-1-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + '.csv', index=False)

sarscov1baseline = 0.1301183

print(sarscov1baseline)


##### This baseline peak-normalized N IgG antibody level for SARS-CoV-1 comes from the ancestral and descendent states analysis that used the baselines for the human "seasonal" coronaviruses to estimate the baselines for the zoonotic coronaviruses.

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

##### This value of lambda goes into the ancestral and descendent states analysis estimating a and b for the zoonotic coronaviruses.

sarscov1halflife = np.log(2) / 0.0015071475711347397 

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

lambdaForPlot = 0.00150715
result = 0
plotBaseLines = list()

for days in range(0, len(sarscov1antibodytimecourse)):
    exponentValue = -lambdaForPlot * days
    result = sarscov1baseline + (1 - sarscov1baseline) * math.exp(exponentValue)
    plotBaseLines.append(result)
    
    
from matplotlib import pyplot as plt

plt.plot(sarscov1antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('SARS-CoV-1 with fixed baseline')

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

f2.savefig("SARS-CoV-1_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

print(len(sarscov1antibodytimecourse))

sarscov1antibodytimecourseplusexp = sarscov1antibodytimecourse
day = len(sarscov1antibodytimecourse)

print(day)

print(len(sarscov1antibodytimecourseplusexp))

#while day < 4393:
#Small file 86 records
#while day < 89:
maxDays = 2
while day < maxDays:
    day = day + 1
    print("day : " + str(day))
    exponentValueBaseLine = -0.0017578418303445613
    #tempValue = sarscov1baseline + (sarscov1antibodytimecourseplusexp[day-1] - sarscov1baseline) * math.exp(exponentValueBaseLine)
    tempValue = sarscov1baseline + (sarscov1antibodytimecourseplusexp[0] - sarscov1baseline) * math.exp(exponentValueBaseLine)
    print('tempValue' + str(tempValue))
    sarscov1antibodytimecourseplusexp.append(tempValue)
    
    
print(len(sarscov1antibodytimecourseplusexp))

import matplotlib.pyplot as plt2
plt2.plot(sarscov1antibodytimecourseplusexp)
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

##### SARS-CoV-1 Probability of Infection | Antibody OD

##### These values for a and b come from the ancestral and descendent states analysis to determine them for the zoonotic coronaviruses given baselines and declines for all viruses, and a's and b's for the "seasonal" coronaviruses.

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


##### The alternate values for a and b above come from our results using different approaches toward building the molecular evolutionary tree of the coronaviruses and toward building the time tree of the coronaviruses (see Supplement).

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

sarscov1probinfList1 = populatesarscov1probinfList(-4.476207, -16.10749)

plt.plot(sarscov1probinfList1, color='red')
plt.ylim(0,0.013)
plt.xlim(0.0, 1.0)
plt.xlabel("Peak-normalized IgG antibody OD")

sarscov1probinfList2 = populatesarscov1probinfList(-4.654767, -14.63486)
sarscov1probinfList3 = populatesarscov1probinfList(-4.102592, -18.02885)
sarscov1probinfList4 = populatesarscov1probinfList(-4.046638, -18.74478)
sarscov1probinfList5 = populatesarscov1probinfList(-4.298427, -16.88932)

plt.plot(sarscov1probinfList2, color = 'green')
plt.plot(sarscov1probinfList3, color = 'lightblue')
plt.plot(sarscov1probinfList4, color = 'purple')
plt.plot(sarscov1probinfList5, color = 'blue')

plt.ylabel("Daily probability of infection")
plt.title('SARS-CoV-1')
#plt.show()

plt.title('SARS-CoV-1')
plt.show()

##### The alternate values for a and b above come from our results using different approaches toward building 
##### the molecular evolutionary tree of the coronaviruses and toward building the time tree of 
##### the coronaviruses (see Supplement).

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

f5.savefig("SARS-CoV-1_PrInfs-by-nOD.pdf", bbox_inches='tight')

##### SARS-CoV-1 Probability of No Reinfection Time Course

sarscov1probnoreinfectiontimecourse = list()
day = 0
sarscov1probnoreinfectiontimecourse.append(1.0)

#Not working
#sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1-sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day+1])))
sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1-sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day])))

print(sarscov1probnoreinfectiontimecourse)

print(len(sarscov1antibodytimecourse))

#### Need to change days to 4393
while day < 1:
    
    if(day < len(sarscov1antibodytimecourse)):
        #sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1 - sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day + 1])))
        sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1 - sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[0])))
    else:
        #sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1 - sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[1392])))
        sarscov1probnoreinfectiontimecourse.append(sarscov1probnoreinfectiontimecourse[day] * (1 - sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[0])))
            
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

len(sarscov1antibodytimecourseplusexp)

len(sarscov1probnoreinfectiontimecourse)

sarscov1probreinfection = list()

day = 0

#while day < 4393:
while day < 1:
    day = day + 1  
        
    #if day < 1392:
    if day < 1:
        tempValue = sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[day]) * sarscov1probnoreinfectiontimecourse[day]
    else:
        tempValue = sarscov1probinfgivenaod(sarscov1antibodytimecourseplusexp[0])* sarscov1probnoreinfectiontimecourse[0]
    
    sarscov1probreinfection.append(tempValue)      
    
import matplotlib.pyplot as plt7
plt7.plot(sarscov1probreinfection)
plt7.ylim(0, 0.0002)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('SARS-CoV-1') 

import matplotlib.pyplot as plt7

f = plt7.figure()

plt7.plot(sarscov1probreinfection)
plt7.ylim(0, 0.0002)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('SARS-CoV-1') 
f.savefig("SARS-CoV-1_PrInf-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-2

##### Input Data: Average peak-normalized ELIZA ODs for SARS-CoV-2 N IgG antibodies converted to IgG antibodies from from Gudbjartsson et al. 2020, "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine. using a linear regression model of N and IgG antibody data from Li et al. 2006, "Long-Term Persistence of Robust Antibody and Cytotoxic T Cell Responses in Recovered Patients Infected with SARS Coronavirus", PLoS ONE. 

##### days : 35, 48, 70, 94, 109. OD : 0.3500076, 0.3397884, 0.3380852, 0.29806, 0.302318

sarscov2data = [[1, 35, "SARSCoV2", 0.3500076/0.3500076, True, "NA", "NA"],
                [2, 48, "SARSCoV2", 0.3397884/0.3500076, False, (48 - 35), (0.3397884 - 0.3500076 / 0.3500076)/(48 - 35)],
                [3, 70, "SARSCoV2", 0.3380852/0.3500076, False, (70 - 48), (0.3380852 - 0.3397884 / 0.3500076)/(70 - 48)], 
                [4, 94, "SARSCoV2", 0.29806/0.3500076, False, (95 - 70), (0.29806 - 0.3380852 / 0.3500076)/(94 - 70)]]      


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
            valueList.append(sarscov2data[index-1][6])
        #else:
         #   valueList.append("None")
        if (len(valueList) != 0 ):
            sarscov2paddedmeanwaning["{0:.3f}".format(aod)] = valueList            
        index = index + 1    
    aod = aod + 0.02

print(len(sarscov2paddedmeanwaning))

print(sarscov2paddedmeanwaning)

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

##### This baseline peak-normalized N IgG antibody level for SARS-CoV-2 comes from the ancestral and descendent states analysis that used the baselines for the human "seasonal" coronaviruses to estimate the baselines for the zoonotic coronaviruses.

sarscov2baseline = 0.1301183

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

##### This value of lambda goes into the ancestral and descendent states analysis estimating a and b for the zoonotic coronaviruses.

sarscov2halflife = np.log(2) / 0.003746311184177791 

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
maxDays = 6
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

print(len(sarscov2antibodytimecourseplusexp))

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

f.savefig("SARS-CoV-2_nOD-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-2 Probability of Infection | Antibody OD

##### These values for a and b come from the ancestral and descendent states analysis to determine them for the zoonotic coronaviruses given baselines and declines for all viruses, and a's and b's for the "seasonal" coronaviruses.

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

##### The alternate values for a and b above come from our results using different approaches toward building the molecular evolutionary tree of the coronaviruses and toward building the time tree of the coronaviruses (see Supplement).

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

#### Need to change days to 4393
while day < 2:
    
    if(day < len(sarscov2antibodytimecourse)):
        sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day + 1])))
    else:
        #sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1392])))
        sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1])))
            
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

f15.savefig("SARS-CoV-2_PnorInfTimecourse-by-" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')


sarscov2probreinfection = list()
day = 0;

#while day < 4393:
while day < 85:
    day = day + 1  
        
    #if day < 1392:
    if day < 2:
        tempValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day]) * sarscov2probnoreinfectiontimecourse[day]
    else:       
        #empValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1392])  * sarscov2probnoreinfectiontimecourse[day]
        tempValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1])  * sarscov2probnoreinfectiontimecourse[1]
    
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