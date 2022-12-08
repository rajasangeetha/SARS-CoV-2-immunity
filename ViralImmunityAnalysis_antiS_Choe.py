##### Durability of Immunity against reinfection of SARS-CoV-2

##### Author: Sangeetha Vijayam
##### Reviwer: Jeffrey P. Townsend, Alex Dornburg, Hayley B. Hassler
##### Date: 10/11/2022

####### Edridge et al. full dataset

####### Input Data: updated_interp_OD _filev7.csv contains processed N IgG antibodies converted to S IgG antibodies from the supplementary dataset of Edridge et al. 2020, "Seasonal coronavirus protective immunity is short-lasting", Nature Medicine. using a linear regression model of N and S IgG antibody data from Gudbjartsson et al. 2020, "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine.

##### 1. Import and get header data

import pandas as pd
import numpy
from datetime import datetime
import math

### Need to change the file 
inputFile = "updated_OD_filev7_org.csv"

df2 = pd.read_csv(inputFile)

print(df2)

headers = list(df2.head(0))
print(headers)

##### Take dataset without headers
noheaderdataset0 = pd.read_csv(inputFile).iloc[0:]
print(noheaderdataset0)

##### Remove Rowws only with NaN => "NA:"
noheaderdataset = noheaderdataset0.dropna(axis = 0, how ='all')
print(noheaderdataset)

print(len(noheaderdataset))

fulldataset = numpy.transpose(noheaderdataset) 
fulldataset = numpy.transpose(fulldataset) 
print(fulldataset)


##### OC43
##### Input Data: ELIZA ODs for HCoV-OC43 N IgG antibodies converted to S IgG antibodies from the supplementary 
##### dataset of Edridge et al. 2020, "Seasonal coronavirus protective immunity is short-lasting", Nature Medicine. 
##### using a linear regression model of N and S IgG antibody data from Gudbjartsson et al. 2020, 
##### "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine.

df = pd.DataFrame(fulldataset)
print(df)
oc43data = df[df.isin(["OC43"]).any(axis=1)]
print(oc43data)

oc43peakdata = df2[df2.isin(["TRUE"]).any(axis=1)]
print(oc43peakdata)

lloc43length = len(oc43data)
print(lloc43length)

##### OC43 Probability of Infection | Antibody OD
######## ToDo: Need to do Full Simplify and PowerExpand
##### OC43 Waning of Antibody OD
print(oc43data)

oc43dataArr = oc43data.to_numpy()

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
oc43paddedmeanwaning = dict();

index = 0

while aod <= 2.9:       
    index  = 2
    valueList = list()
    
    while index <= lloc43length:    
        #print("aod :" + str(aod))
        #print("index : " + str(index))
        print(oc43dataArr[index-2])
        
        if ((oc43dataArr[index-2][9] >= aod and aod >= oc43dataArr[index-2][9]) or 
            (oc43dataArr[index-2][9] <= aod and aod <= oc43dataArr[index-2][9]) ) :
            if (oc43dataArr[[index-1, 10]] != "NA" and re.match( r'^(F).*(E)$', oc43dataArr[index-2][6] )):
                if oc43dataArr[index-1][8].isnumeric() == True: 
                    valueList.append(oc43dataArr[index-1][8])
                else:
                    valueList.append(0)
        average = grade_avg(valueList)
        print("average : " + str(average))
        oc43paddedmeanwaning["{0:.3f}".format(aod)] = average
        
        index = index + 1    
    aod = aod + 0.05
    
 print(oc43paddedmeanwaning)
 
for key in dict(oc43paddedmeanwaning):
    print(key)

(k := next(iter(oc43paddedmeanwaning)), oc43paddedmeanwaning.pop(k))
(k := next(iter(oc43paddedmeanwaning)), oc43paddedmeanwaning.pop(k))


updict = {'0.05' : '0', '0': '0'}

# ** operator for packing and unpacking items in order
oc43paddedmeanwaningNew = {**updict, **oc43paddedmeanwaning}
print(oc43paddedmeanwaningNew)

import numpy as np

aodList = list(oc43paddedmeanwaningNew.keys())
print(aodList)

arraod = np.array(aodList)
arraod = arraod.astype(float)

covdataList = list(oc43paddedmeanwaningNew.values())
print(covdataList)

arrcov = np.array(covdataList).squeeze()

from scipy.interpolate import interp1d
print(len(aodList))
print(len(covdataList))

y_interpolation = interp1d(arraod, arrcov)

interpolate_x_new = 1.0


print("Value of Y at x = {} is".format(interpolate_x_new), y_interpolation(interpolate_x_new))
      
      
oc43antibodytimecourse = list()

oc43antibodytimecourse.append(1.0)
 
print(str(oc43antibodytimecourse))

##### Populating a dictionary for the interpolate values

day = 0
print("Populating a dictionary for the interpolate values")

while oc43antibodytimecourse[day] >= 0.1:
    #print("For day : " +  str(day))
    #print("oc43antibodytimecourse[day] : " + str(oc43antibodytimecourse[day]))
    #print("oc43antibodytimecourse : " + str(oc43antibodytimecourse))
    #print("inerpolate of sarscov1paddedmeanwaning")
    #print(y_interpolation(oc43antibodytimecourse[day]))
    #print("Added value")
    oc43antibodytimecourse.append(oc43antibodytimecourse[day] + y_interpolation(oc43antibodytimecourse[day]))    
    day = day + 1
    
print(oc43antibodytimecourse)

from matplotlib import pyplot as plt
plt.plot(oc43antibodytimecourse)
plt.ylim(0,1.0)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-OC43')
plt.show()


import matplotlib.pyplot as plt
import math

f = plt.figure()
plt.plot(oc43antibodytimecourse)
plt.ylim(0,1.0)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-OC43')
plt.show()

f.savefig("HCoV-OC43_AntibodyTimecourse" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

###### Export to Excel

from datetime import datetime

df3 = pd.DataFrame(oc43antibodytimecourse, columns=["Antibody Time Course"])
df3.to_csv('OC43-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + '.csv', index=False)

len(oc43antibodytimecourse)

import math

def calculatelambda(sarscov1baseline, lambdaValue, index):
    exponentValue = -lambdaValue * index
    result = sarscov1baseline + (1 - sarscov1baseline) * math.exp(exponentValue)
    return result
    
oc43lsfunc = 0
baseline = 0.1
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(oc43antibodytimecourse))):
    print(oc43antibodytimecourse[index])
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.005
    
    ##Added 11/28
    oc43lsfunc = (oc43antibodytimecourse[index]) - math.pow((calculatelambda(baseline, lambdaValue, index)), 2)     
    print(oc43lsfunc)
    
    if oc43lsfunc < currentMinValue:
        currentMinValue = oc43lsfunc      

print(currentMinValue)

oc43halflife = np.log(2) / 0.0042019881376685235

print(oc43halflife)

##### Added as part of python code

plt.plot(oc43antibodytimecourse, color='red')
plt.ylim(0,1.0)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-OC43 with fixed baseline')
plt.show()

##### Added as part of python code
f2 = plt.figure()
plt.plot(oc43antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-OC43 with fixed baseline')
plt.show()

f2.savefig("HCoV-OC43_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lambdaForPlot = 0.00420101
baseline = 0.0996472
result = 0
plotBaseLines = list()

for days in range(0, len(oc43antibodytimecourse)):
    exponentValue = -lambdaForPlot * days
    result = baseline + (1 - baseline) * math.exp(exponentValue)
    plotBaseLines.append(result)
    
print(plotBaseLines)


from matplotlib import pyplot as plt

plt.plot(oc43antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-OC43 with fixed baseline')
#plt.show()

plt.title('HCoV-OC43')
plt.show()


##### Added as part of python code

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

len(oc43antibodytimecourse)

oc43antibodytimecourseplusexp = oc43antibodytimecourse
day = len(oc43antibodytimecourse)

print(day)

#while day < 4393 - Length[oc43antibodytimecourse]:
#Small file 86 records
#while day < 89:
maxDays = 2
while day < maxDays:
    day = day + 1
    print("day : " + str(day))
    exponentValueBaseLine = -0.00420101
    tempValue = 0.09964716261621707 + (oc43antibodytimecourseplusexp[day-1] - 0.09964716261621707) * math.exp(exponentValueBaseLine)
    print('tempValue' + str(tempValue))
    oc43antibodytimecourseplusexp.append(tempValue)
    
print(oc43antibodytimecourseplusexp)

import matplotlib.pyplot as plt2
plt2.plot(oc43antibodytimecourseplusexp)
plt2.ylim(0,1.1)
plt2.xlabel("Day")
plt2.ylabel("Peak-normalized antibody OD")
plt2.title('HCoV-OC43')

import matplotlib.pyplot as plt3

f = plt3.figure()
plt3.plot(oc43antibodytimecourseplusexp)

plt3.ylim(0,1.1)
plt3.xlabel("Day")
plt3.ylabel("Peak-normalized antibody OD")
plt3.title('HCoV-OC43')
plt3.show()

f.savefig("HCoV-OC43_nOD-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lengthabtc = len(oc43antibodytimecourse)
print(lengthabtc)

##### HCoV-OC43 Probability of No Reinfection Time Course

######## ToDo: Need to update based on Full Simplify value
a = 1
b = 1

import math

def oc43probinfgivenaod(aod, a, b):
    exponentValue = -(a + b * aod)   
    result = 1 / (1 + math.exp(exponentValue))
    return result

oc43probnoreinfectiontimecourse = list()
day = 0
oc43probnoreinfectiontimecourse.append(1.0)

oc43probnoreinfectiontimecourse.append(oc43probnoreinfectiontimecourse[day] * (1-oc43probinfgivenaod(oc43antibodytimecourse[day+1], a, b)))

print(oc43probnoreinfectiontimecourse)

#### Need to change days to 4393
while day < lengthabtc:  
    print("day : " +str(day))
    if(day < lengthabtc):
        #Not working
        #oc43probnoreinfectiontimecourse.append(oc43probnoreinfectiontimecourse[day] * (1 - oc43probinfgivenaod(oc43antibodytimecourse[day + 1], a, b)))
        oc43probnoreinfectiontimecourse.append(oc43probnoreinfectiontimecourse[day] * (1 - oc43probinfgivenaod(oc43antibodytimecourse[day], a, b)))
    else:
        oc43probnoreinfectiontimecourse.append(oc43probnoreinfectiontimecourse[day] * (1 - oc43probinfgivenaod(oc43antibodytimecourse[lengthabtc], a, b)))
            
    day = day + 1   
    
print(oc43probnoreinfectiontimecourse)

plt.plot(oc43probnoreinfectiontimecourse, color='blue')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Probability of no reinfection")
plt.title('HCoV-OC43')
plt.show()

import matplotlib.pyplot as plt6

f6 = plt6.figure()
plt6.plot(oc43probnoreinfectiontimecourse, color = 'blue')

#plt6.ylim(0,0.007)
#plt6.xlim(0.0, 1.0)

plt6.xlabel("Day")
plt6.ylabel("Probability of no reinfection")
plt6.title('HCoV-OC43')
plt6.show()

f6.savefig("HCoV-NL63_PnorInfTimecourse-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

oc43probreinfection = list()

day = 0

#while day < 4393:
while day < 2:
    day = day + 1  
        
    #if day < 1392:
    if day < 2:
        tempValue = oc43probinfgivenaod(oc43antibodytimecourse[day], a, b) * oc43probnoreinfectiontimecourse[day]
    else:
        #tempValue = oc43probinfgivenaod(oc43antibodytimecourse[lengthabtc], a, b)* oc43probnoreinfectiontimecourse[day]
        tempValue = oc43probinfgivenaod(oc43antibodytimecourse[1], a, b)* oc43probnoreinfectiontimecourse[day]
    
    oc43probreinfection.append(tempValue)      
    
import matplotlib.pyplot as plt7
plt7.plot(oc43probreinfection)
plt7.ylim(0, 1.0)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('HCoV-OC43') 

import matplotlib.pyplot as plt7

f = plt7.figure()

plt7.plot(oc43probreinfection)
plt7.ylim(0, 1.0)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('HCoV-OC63') 
f.savefig("HCoV-OC43_PrInf-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### NL63 Waning of Antibody OD
###### Input Data: ELIZA ODs for HCoV-NL63 N IgG antibodies converted to S IgG antibodies from the 
###### supplementary dataset of Edridge et al. 2020, "Seasonal coronavirus protective immunity is short-lasting", 
###### Nature Medicine. using a linear regression model of N and S IgG antibody data from Gudbjartsson et al. 2020, 
###### "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine.

nl63data = df[df2.isin(["NL63"]).any(axis=1)]
print(nl63data)

print(len(nl63data))

####### To Do: need to calculate llnl63lrdaily to get FullSimplify
a = 1
b = 1

import math

def nl63probinfgivenaod(aod, a, b):
    exponentValue = -(a + b * aod)   
    result = 1 / (1 + math.exp(exponentValue))
    return result  
    
    
nl63dataArr = nl63data.to_numpy()

import re

aod = 0
nl63paddedmeanwaning = dict();

index = 0

while aod <= 2.9:       
    index  = 2
    valueList = list()
    
    while index <= lloc43length:    
        #print("aod :" + str(aod))
        #print("index : " + str(index))
        print(nl63dataArr[index-2])
        
        if ((nl63dataArr[index-2][9] >= aod and aod >= nl63dataArr[index-2][9]) or 
            (nl63dataArr[index-2][9] <= aod and aod <= nl63dataArr[index-2][9]) ) :
            if (nl63dataArr[[index-1, 10]] != "NA" and re.match( r'^(F).*(E)$', nl63dataArr[index-2][6] )):
                if nl63dataArr[index-1][8].isnumeric() == True: 
                    valueList.append(oc43dataArr[index-1][8])
                else:
                    valueList.append(0)
        average = grade_avg(valueList)
        print("average : " + str(average))
        nl63paddedmeanwaning["{0:.3f}".format(aod)] = average
        
        index = index + 1    
    aod = aod + 0.05

print(nl63paddedmeanwaning)

for i in range(0, 5):
    (k := next(iter(nl63paddedmeanwaning)), nl63paddedmeanwaning.pop(k))

print(nl63paddedmeanwaning)

updict = {'0.20':'0', '0.15':'0', '0.10':'0', '0.05':'0', '0':'0'}

# ** operator for packing and unpacking items in order
nl63paddedmeanwaningNew = {**updict, **nl63paddedmeanwaning}
print(nl63paddedmeanwaningNew)

import numpy as np

aodList = list(nl63paddedmeanwaningNew.keys())
print(aodList)

arraod = np.array(aodList)
arraod = arraod.astype(float)

covdataList = list(nl63paddedmeanwaningNew.values())
print(covdataList)

arrcov = np.array(covdataList).squeeze()

from scipy.interpolate import interp1d
print(len(aodList))
print(len(covdataList))

y_interpolation = interp1d(arraod, arrcov)

interpolate_x_new = 1.0


print("Value of Y at x = {} is".format(interpolate_x_new),
      y_interpolation(interpolate_x_new))
      
nl63antibodytimecourse = list()

nl63antibodytimecourse.append(1.0)

print(str(nl63antibodytimecourse))

##### Populating a dictionary for the interpolate values


day = 0
print("Populating a dictionary for the interpolate values")

while nl63antibodytimecourse[day] >= 0:
    print("For day : " +  str(day))
    print("nl63antibodytimecourse[day] : " + str(nl63antibodytimecourse[day]))
    print("nl63antibodytimecourse : " + str(nl63antibodytimecourse))
    print("inerpolate of nl63paddedmeanwaningNew")
    print(y_interpolation(nl63antibodytimecourse[day]))
    print("Added value")
    nl63antibodytimecourse.append(nl63antibodytimecourse[day] + y_interpolation(nl63antibodytimecourse[day]))    
    day = day + 1
    
print(nl63antibodytimecourse)

from matplotlib import pyplot as plt
plt.plot(nl63antibodytimecourse)
plt.ylim(0,1.0)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-NL63')
plt.show()

import matplotlib.pyplot as plt
import math

f = plt.figure()
plt.plot(nl63antibodytimecourse)
plt.ylim(0,1.0)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-NL63')
plt.show()

f.savefig("HCoV-NL63_AntibodyTimecourse" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

nl63lsfunc = 0
baseline = 0.1
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(nl63antibodytimecourse))):
    print(nl63antibodytimecourse[index])
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.005
    
    ##Added 11/28
    nl63lsfunc = (nl63antibodytimecourse[index]) - math.pow((calculatelambda(baseline, lambdaValue, index)), 2)     
    print(oc43lsfunc)
    
    if nl63lsfunc < currentMinValue:
        currentMinValue = nl63lsfunc    
        
print(currentMinValue)

nl63halflife  = np.log(2) / 0.0023166874396891007

print(nl63halflife )

plt.plot(nl63antibodytimecourse, color='red')
plt.ylim(0,1.0)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-NL63 with fixed baseline')
plt.show()

f2 = plt.figure()
plt.plot(nl63antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-NL63 with fixed baseline')
plt.show()

f2.savefig("HCoV-NL63_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

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

plt.plot(nl63antibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-NL63 with fixed baseline')
#plt.show()

plt.title('HCoV-NL63')
plt.show()

f2 = plt.figure()
plt.plot(plotBaseLines)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-NL63')
plt.show()

f2.savefig("HCoV-NL63_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

nl63probnoreinfectiontimecourse = list()
day = 0
nl63probnoreinfectiontimecourse.append(1.0)

nl63probnoreinfectiontimecourse.append(nl63probnoreinfectiontimecourse[day] * (1-nl63probinfgivenaod(nl63antibodytimecourse[day+1], a, b)))

print(nl63probnoreinfectiontimecourse)

#### Need to change days to 4393
while day < 85:
    
    if(day < len(nl63probnoreinfectiontimecourse)):        
        nl63probnoreinfectiontimecourse.append(nl63probnoreinfectiontimecourse[day] * (1 - nl63probinfgivenaod(nl63probnoreinfectiontimecourse[day + 1], a, b)))        
    else:        
        nl63probnoreinfectiontimecourse.append(nl63probnoreinfectiontimecourse[day] * (1 - nl63probinfgivenaod(nl63probnoreinfectiontimecourse[1392], a, b)))
            
    day = day + 1   
    
print(nl63probnoreinfectiontimecourse)

plt.plot(nl63probnoreinfectiontimecourse, color='blue')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Probability of no reinfection")
plt.title('HCoV-NL63')
plt.show()

import matplotlib.pyplot as plt6

f6 = plt6.figure()
plt6.plot(nl63probnoreinfectiontimecourse, color = 'blue')

#plt6.ylim(0,0.007)
#plt6.xlim(0.0, 1.0)

plt6.xlabel("Day")
plt6.ylabel("Probability of no reinfection")
plt6.title('HCoV-NL63')
plt6.show()

f6.savefig("HCoV-NL63_PnorInfTimecourse-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

nl63probreinfection = list()

day = 0

#while day < 4393:
while day < 4393:
    day = day + 1  
        
    #if day < 1392:
    if day < 2:
        tempValue = nl63probinfgivenaod(oc43antibodytimecourse[day], a, b) * nl63probnoreinfectiontimecourse[day]
    else:
        #Not working
        #tempValue = oc43probinfgivenaod(oc43antibodytimecourse[lengthabtc], a, b)* oc43probnoreinfectiontimecourse[day]
        tempValue = nl63probinfgivenaod(oc43antibodytimecourse[1], a, b)* nl63probnoreinfectiontimecourse[1]
    
    nl63probreinfection.append(tempValue)      
    
    
print(nl63probreinfection)

import matplotlib.pyplot as plt7
plt7.plot(nl63probreinfection)
plt7.ylim(0, 1.0)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('HCoV-NL63') 

import matplotlib.pyplot as plt7

f = plt7.figure()

plt7.plot(nl63probreinfection)
plt7.ylim(0, 1.0)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('HCoV-OC43') 
f.savefig("HCoV-NL63_PrInf-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### 229E Waning of Antibody OD
####### Input Data: ELIZA ODs for HCoV-229E N IgG antibodies converted to S IgG antibodies from the supplementary dataset 
####### of Edridge et al. 2020, "Seasonal coronavirus protective immunity is short-lasting", Nature Medicine. 
####### using a linear regression model of N and S IgG antibody data from Gudbjartsson et al. 2020, 
####### "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine.

x229edata = df2[df2.isin(["229E"]).any(axis=1)]

print(x229edata)

print(x229edata.iloc[:1])

llx229elength = len(x229edata)

print(llx229elength)

##### Probability of Infection | Antibody OD

##### Probability of Infection | Antibody OD

###### llx229elrdaily = yet to do

##### 229E Waning of Antibody OD

print(x229edata)

x229edataArr = x229edata.to_numpy()

import re

aod = 0
x229epaddedmeanwaning = dict();

index = 0

while aod <= 1.5:       
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
    
for key in dict(x229epaddedmeanwaning):
    print(key)
    
(k := next(iter(x229epaddedmeanwaning)), x229epaddedmeanwaning.pop(k))

updict = {'0.05':'0', '0': '0'}

# ** operator for packing and unpacking items in order
x229epaddedmeanwaningNew = {**updict, **x229epaddedmeanwaning}
print(x229epaddedmeanwaningNew)

x229emeanwaning = x229epaddedmeanwaningNew

import numpy as np

aodList = list(x229emeanwaning.keys())
print(aodList)

arraod = np.array(aodList)
arraod = arraod.astype(float)

covdataList = list(x229emeanwaning.values())
print(covdataList)

arrcov = np.array(covdataList).squeeze()

from scipy.interpolate import interp1d
print(len(aodList))
print(len(covdataList))

y_interpolation = interp1d(arraod, arrcov)

x229eantibodytimecourse = list()

x229eantibodytimecourse.append(1.0)

print(str(x229eantibodytimecourse))

##### Populating a dictionary for the interpolate values

day = 0
print("Populating a dictionary for the interpolate values")

while x229eantibodytimecourse[day] >= 0.1:
    print("For day : " +  str(day))
    print("x229eantibodytimecourse[day] : " + str(x229eantibodytimecourse[day]))
    print("x229eantibodytimecourse : " + str(x229eantibodytimecourse))
    print("inerpolate of x229emeanwaning")
    print(y_interpolation(x229eantibodytimecourse[day]))
    print("Added value")
    x229eantibodytimecourse.append(x229eantibodytimecourse[day] + y_interpolation(x229eantibodytimecourse[day]))    
    day = day + 1
    
print(x229eantibodytimecourse)

from matplotlib import pyplot as plt
plt.plot(x229eantibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-229E')
plt.show()

import matplotlib.pyplot as plt
import math

f = plt.figure()
plt.plot(x229eantibodytimecourse)
plt.ylim(0,1.0)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('HCoV-229E')
plt.show()

f.savefig("HCoV-229E_AntibodyTimecourse" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

###### Export to Excel

from datetime import datetime

df3 = pd.DataFrame(x229eantibodytimecourse, columns=["Antibody Time Course"])
df3.to_csv('HCoV-229E-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + '.csv', index=False)

print(len(x229eantibodytimecourse))

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
    x229elsfunc = (x229eantibodytimecourse[index]) - math.pow((calculatelambda(baseline, lambdaValue, index)), 2)     
    print(x229elsfunc)
    
    if x229elsfunc < currentMinValue:
        currentMinValue = x229elsfunc      
        
print(currentMinValue)

x229ehalflife = np.log(2) / 0.005241500716948157

print(x229ehalflife)

##### Added as part of python code

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.0)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E with fixed baseline')
plt.show()

##### Added as part of python code


f2 = plt.figure()
plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E with fixed baseline')
plt.show()

f2.savefig("HCoV-229E_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

lambdaForPlot = 0.005241500716948157
baseline = 0.0999574
result = 0
plotBaseLines = list()

for days in range(0, len(x229eantibodytimecourse)):
    exponentValue = -lambdaForPlot * days
    result = baseline + (1 - baseline) * math.exp(exponentValue)
    plotBaseLines.append(result)
    
from matplotlib import pyplot as plt

plt.plot(x229eantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('HCoV-229E with fixed baseline')
#plt.show()

plt.title('HCoV-229E')
plt.show()

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
    exponentValueBaseLine = -0.005241500716948157
    tempValue = 0.09995743181844773 + (x229eantibodytimecourseplusexp[day-1] - 0.09995743181844773) * math.exp(exponentValueBaseLine)
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
plt3.title('HCoV-229E')
plt3.show()

f.savefig("HCoV-229E_nOD-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

import math

a = 1
b = 1

def x229eprobinfgivenaod(aod, a, b):
    exponentValue = -(a + b * aod)   
    result = 1 / (1 + math.exp(exponentValue))
    return result
    
x229eprobnoreinfectiontimecourse = list()
day = 0
x229eprobnoreinfectiontimecourse.append(1.0)

x229eprobnoreinfectiontimecourse.append(x229eprobnoreinfectiontimecourse[day] * (1-x229eprobinfgivenaod(x229eantibodytimecourse[day+1], a, b)))

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
        tempValue = x229eprobinfgivenaod(x229eantibodytimecourse[day], a, a) * x229eantibodytimecourse[day]
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


mersdata = [[1, 50, "MERS", 2.67561878/2.67561878, True, "NA", "NA"],
                [2, 75, "MERS", 2.5938248/2.67561878,  False, (75 - 50), ((2.5938248 - 2.67561878)/2.67561878)/(75 - 50)],
                [3, 100, "MERS", 2.49390066/2.67561878, False, (100 - 75),  ((2.49390066 - 2.5938248)/2.67561878)/(100 - 75)],
                [4, 125, "MERS", 2.40265641/2.67561878, False, (125 - 100),  ((2.40265641 - 2.49390066)/2.67561878)/(125 - 100)],
                [5, 150, "MERS", 2.31624222/2.67561878, False, (150 - 125),  ((2.31624222 - 2.40265641)/2.67561878)/(150 - 125)],
                [6, 175, "MERS", 2.22800926/2.67561878, False, (175 - 150),  ((2.22800926 - 2.31624222)/2.67561878)/(175 - 150)],
                [7, 200, "MERS", 2.13682784/2.67561878, False, (200 - 175),  ((2.13682784 - 2.22800926)/2.67561878)/(200 - 175)],
                [8, 225, "MERS", 2.04875422/2.67561878, False, (225 - 200),  ((2.04875422 - 2.13682784)/2.67561878)/(225 - 200)],
                [9, 250, "MERS", 2.00596085/2.67561878, False, (250 - 225),  ((2.00596085 - 2.04875422)/2.67561878)/(250 - 225)],
                [10, 275, "MERS", 1.96931203/2.67561878, False, (275 - 250),  ((1.96931203 - 2.00596085)/2.67561878)/(275 - 250)],
                [11, 300, "MERS", 1.93628719/2.67561878, False, (300 - 275),  ((1.93628719 - 1.96931203)/2.67561878)/(300 - 275)],
                [12, 325, "MERS", 1.88343661/2.67561878, False, (325 - 300),  ((1.88343661 - 1.93628719)/2.67561878)/(325 - 300)]]
                
##### MERS Waning of Antibody OD

mersdatalength = len(mersdata)

print(mersdatalength)

aod3 = 0.725
merspaddedmeanwaning  = dict();

index = 0

while aod3 <= 1:    
    #print("first loop")
   # print("{0:.3f}".format(aod3))
    
    #sarscov1paddedmeanwaning.append("{0:.3f}".format(aod))
    
    index  = 2
    valueList = list()
    #print("index" + str(index))
    
    while index <= mersdatalength:    
        #print(str(index))
        print("aod3:" + str(aod3))

        if ( (mersdata[index-2][3] >= aod3 and aod3 >= mersdata[index-1][3]) or 
             (mersdata[index-2][3] <= aod3 and aod3 <= mersdata[index-1][3]) ) :
            print(mersdata[index-1][6])
            valueList.append(mersdata[index-1][6])

        #else:
         #   valueList.append("None")
        merspaddedmeanwaning["{0:.3f}".format(aod3)] = valueList
            
        index = index + 1
        
    
    #print("second loop")
    aod3 = aod3 + 0.025

print(len(merspaddedmeanwaning))

for key in merspaddedmeanwaning:
    print(str(key) + " " + str(merspaddedmeanwaning[key]))
    
import numpy as np

aodList = list(merspaddedmeanwaning.keys())
aodList

arraod = np.array(aodList)
arraod = arraod.astype(float)

covdataList = list(merspaddedmeanwaning.values())
covdataList

arrcov = np.array(covdataList).squeeze()

from scipy.interpolate import interp1d
print(len(aodList))
print(len(covdataList))

y_interpolation = interp1d(arraod, arrcov)

interpolate_x_new = 0.725


print("Value of Y at x = {} is".format(interpolate_x_new),
      y_interpolation(interpolate_x_new))
      
      
mersantibodytimecourse = list()

mersantibodytimecourse.append(1.0)

print(mersantibodytimecourse)

print(len(mersantibodytimecourse))

print(len(merspaddedmeanwaning))

day = 0
print(mersantibodytimecourse[day])
print(len(mersantibodytimecourse))

while (day < len(merspaddedmeanwaning) and mersantibodytimecourse[day] >= 0.725):  
    print(day)
    print(mersantibodytimecourse[day])    
    #sarscov1antibodytimecourse.append(2.0)
    print(mersantibodytimecourse[day])
    print("merspaddedmeanwaning")
    print(merspaddedmeanwaning)
    print("mersantibodytimecourse[day]")
    print(mersantibodytimecourse[day])
    #print("sarscov1paddedmeanwaning3[sarscov1antibodytimecourse[day]]")
    #sarscov1antibodytimecourse.append(2.0)
    #print(sarscov1paddedmeanwaning3[sarscov1antibodytimecourse[day]])
    print("inerpolate of merspaddedmeanwaning")
    #print(y_interpolation(mersantibodytimecourse[day]))
    print("Added value")
    #test = mersantibodytimecourse[day] + y_interpolation(mersantibodytimecourse[day])
    #print(str())
    #Not working
    #mersantibodytimecourse.append(mersantibodytimecourse[day] + y_interpolation(mersantibodytimecourse[day]))
    mersantibodytimecourse.append(mersantibodytimecourse[day] + 0)
    print("mersantibodytimecourse")
    print(mersantibodytimecourse)

    print("day")
    day = day + 1    
    
    
print(mersantibodytimecourse)

from matplotlib import pyplot as plt
plt.plot(mersantibodytimecourse)
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized antibody OD")
plt.title('MERS')
plt.show()

f2 = plt.figure()
plt.plot(mersantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('MERS')
plt.show()

f2.savefig("MERS_choe_Antibody-Time-Course" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

###### Export to Excel

df3 = pd.DataFrame(mersantibodytimecourse, columns=["Antibody Time Course"])
df3.to_csv('MERS_choe-Antibody-Time-Course.csv', index=False)

mersbaseline = 0.12992406

print(mersbaseline)

##### This baseline peak-normalized S IgG antibody level for MERS comes from the ancestral and descendent states 
##### analysis that used the baselines for the human "seasonal" coronaviruses to estimate the baselines for the 
##### zoonotic coronaviruses.

import math

def calculatelambda(sarscov1baseline, lambdaValue, index):
    exponentValue = -lambdaValue * index
    result = sarscov1baseline + (1 - sarscov1baseline) * math.exp(exponentValue)
    return result
    
merslsfuncwfixedbaseline = 0
currentMinValue = float('inf')
currentMinValue

for index in range(0, (len(mersantibodytimecourse))):
    print(mersantibodytimecourse[index] )
    ##Need to find the step value and maximum value for this lambdaValue
    lambdaValue = 0.002
    
    ##Added 11/28
    merslsfuncwfixedbaseline = (mersantibodytimecourse[index]) - math.pow((calculatelambda(mersbaseline, lambdaValue, index)), 2)     
    print(merslsfuncwfixedbaseline)
    
    if merslsfuncwfixedbaseline < currentMinValue:
        currentMinValue = merslsfuncwfixedbaseline      
        

print(currentMinValue)

#####This value of lambda goes into the ancestral and descendent states analysis estimating a and b for the zoonotic coronaviruses.
mershalflife = np.log(2) / 0.001658632767524445 

print(mershalflife)

##### Added as part of python code

plt.plot(mersantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('MERS with fixed baseline')
plt.show()

##### Added as part of python code


f2 = plt.figure()
plt.plot(mersantibodytimecourse)
plt.ylim(0,1.1)
#####Added x limit
plt.xlim(0,85)
plt.xlabel("Day")
plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('MERS with fixed baseline')
plt.show()

f2.savefig("MERS_withFixedBaseline" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

print(len(mersantibodytimecourse))

lambdaForPlot = 0.001658632767524445
result = 0
plotBaseLines = list()

for days in range(0, len(mersantibodytimecourse)):
    exponentValue = -lambdaForPlot * days
    result = mersbaseline + (1 - mersbaseline) * math.exp(exponentValue)
    plotBaseLines.append(result)
    
print(plotBaseLines)

from matplotlib import pyplot as plt

plt.plot(mersantibodytimecourse, color='red')
plt.ylim(0,1.1)
plt.xlabel("Day")

plt.plot(plotBaseLines, color = 'blue')

plt.ylabel("Peak-normalized IgG antibody OD")
plt.title('MERS with fixed baseline')
#plt.show()
plt.title('MERS with fixed baseline')
plt.show()


print(len(mersantibodytimecourse))

mersantibodytimecourseplusexp = mersantibodytimecourse
day = len(mersantibodytimecourse)

print(day)

print(len(mersantibodytimecourseplusexp))

#while day < 4393:
#Small file 86 records
#while day < 89:
maxDays = 2
while day < maxDays:
    day = day + 1
    print("day : " + str(day))
    exponentValueBaseLine = -0.001658632767524445
    tempValue = mersbaseline + (mersantibodytimecourseplusexp[day-1] - mersbaseline) * math.exp(exponentValueBaseLine)
    print('tempValue' + str(tempValue))
    mersantibodytimecourseplusexp.append(tempValue)
    
print(mersantibodytimecourseplusexp)

import matplotlib.pyplot as plt2
plt2.plot(mersantibodytimecourseplusexp)
plt2.ylim(0,1.1)
plt2.xlabel("Day")
plt2.ylabel("Peak-normalized antibody OD")
plt2.title('MERS')


import matplotlib.pyplot as plt3

f = plt3.figure()
plt3.plot(mersantibodytimecourseplusexp)

plt3.ylim(0,1.1)
plt3.xlabel("Day")
plt3.ylabel("Peak-normalized antibody OD")
plt3.title('MERS')
plt3.show()

f.savefig("MERS_choe_nOD-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### MERS Probability of Infection | Antibody OD

import math

a = 1
b = 1

def mersprobinfgivenaod(aod, a, b):
    exponentValue = -(a + b * aod)   
    result = 1 / (1 + math.exp(exponentValue))
    return result
    
##### These values for a and b come from the ancestral and descendent states analysis to determine them for the 
##### zoonotic coronaviruses given baselines and declines for all viruses, 
##### and a's and b's for the "seasonal" coronaviruses.

mersprobinfgivenaodList = list()

currentValue = 0;

while( currentValue <= 1):
    currentValue = currentValue + 0.00625
    mersprobinfgivenaodList.append(mersprobinfgivenaod(currentValue, a, b))

print(mersprobinfgivenaodList)

plt.plot(mersprobinfgivenaodList)
plt.ylim(0, 0.007)
plt.xlim(0.0, 1.0)
plt.xlabel("Peak - normalized IgG antibody OD")
plt.ylabel("Daily probability of infection")
plt.title('MERS')
plt.show()

import matplotlib.pyplot as plt4

f2 = plt4.figure()
plt4.plot(mersprobinfgivenaodList)

plt4.ylim(0, 0.007)
plt4.xlim(0.0, 1.0)
plt4.xlabel("Peak - normalized IgG antibody OD")
plt4.ylabel("Daily probability of infection")
plt4.title('MERS')
plt4.show()

f2.savefig("MERS_choe_PrInf-by-nOD_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

import math

def mersprobinfgivenall(a, b, aod):
    exponentValue = (-a) + (-b) * aod
    result = 1 / (1 + math.exp(exponentValue))
    return result

def populatemersprobinfList(a, b):
    mersprobinfList = list()
    
    x = 0;

    while( x <= 1):
        x = x + 0.05
        mersprobinfList.append(mersprobinfgivenall(a, b, x))
    
    return mersprobinfList  
    
from matplotlib import pyplot as plt

mersprobinfList1 = populatemersprobinfList(-4.752999, -12.98306)

plt.plot(mersprobinfList1, color='red')
plt.ylim(0,0.009)
plt.xlim(0.0, 1.0)
plt.xlabel("Peak-normalized IgG antibody OD")

mersprobinfList2 = populatemersprobinfList(-4.1612, -15.91843)
mersprobinfList3 = populatemersprobinfList(-5.563915, -7.994218)
mersprobinfList4 = populatemersprobinfList(-5.908631, -6.627999)
mersprobinfList5 = populatemersprobinfList(-4.829109, -12.16748)

plt.plot(mersprobinfList2, color = 'green')
plt.plot(mersprobinfList3, color = 'lightblue')
plt.plot(mersprobinfList4, color = 'purple')
plt.plot(mersprobinfList5, color = 'blue')

plt.ylabel("Daily probability of infection")
plt.title('MERS')
plt.show()

##### The alternate values for a and b above come from our results using different approaches toward 
##### building the molecular evolutionary tree of the coronaviruses and toward building the time tree of the 
##### coronaviruses (see Supplement).

import matplotlib.pyplot as plt5

f5 = plt5.figure()
plt5.plot(mersprobinfList1, color='red')
plt5.plot(mersprobinfList2, color = 'green')
plt5.plot(mersprobinfList3, color = 'lightblue')
plt5.plot(mersprobinfList4, color = 'purple')
plt5.plot(mersprobinfList5, color = 'blue')

plt5.ylim(0,0.007)
plt5.xlim(0.0, 1.0)

plt5.xlabel("Peak - normalized IgG antibody OD")
plt5.ylabel("Daily probability of infection")
plt5.title('SARS-CoV-1')
plt5.show()

f5.savefig("SARS-CoV-1_PrInfs-by-nOD.pdf", bbox_inches='tight')

##### MERS Probability of No Reinfection Time Course

mersprobnoreinfectiontimecourse = list()
day = 0
mersprobnoreinfectiontimecourse.append(1.0)

mersprobnoreinfectiontimecourse.append(mersprobnoreinfectiontimecourse[day] * (1- mersprobinfgivenaod(mersantibodytimecourseplusexp[day+1], a, b)))

print(mersprobnoreinfectiontimecourse)

print(len(mersantibodytimecourse))

#### Need to change days to 4393
while day < 2:
    
    if(day < len(mersantibodytimecourse)):
        #mersprobnoreinfectiontimecourse.append(mersprobnoreinfectiontimecourse[day] * (1 - mersprobinfgivenaod(mersantibodytimecourseplusexp[day + 1])))
        mersprobnoreinfectiontimecourse.append(mersprobnoreinfectiontimecourse[day] * (1 - mersprobinfgivenaod(mersantibodytimecourseplusexp[1], a, b)))
    else:
        mersprobnoreinfectiontimecourse.append(mersprobnoreinfectiontimecourse[day] * (1 - mersprobinfgivenaod(mersantibodytimecourseplusexp[1], a, b)))
            
    day = day + 1   
    
print(mersprobnoreinfectiontimecourse)

plt.plot(mersprobnoreinfectiontimecourse, color='blue')
#Limit not working
#plt.ylim(0,1.1)
plt.xlabel("Day")
plt.ylabel("Probability of no reinfection")
plt.title('MERS')
plt.show()

import matplotlib.pyplot as plt6

f6 = plt6.figure()
plt6.plot(mersprobnoreinfectiontimecourse, color = 'blue')

#plt6.ylim(0,0.007)
#plt6.xlim(0.0, 1.0)

plt6.xlabel("Day")
plt6.ylabel("Probability of no reinfection")
plt6.title('MERS')
plt6.show()

f6.savefig("MERS_PnorInfTimecourse-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

mersprobreinfection = list()

day = 0

#while day < 4393:
while day < 2:
    day = day + 1  
        
    #if day < 1392:
    if day < 2:
        tempValue = mersprobinfgivenaod(mersantibodytimecourseplusexp[day], a, b) * mersprobnoreinfectiontimecourse[day]
    else:           
        tempValue = mersprobinfgivenaod(mersantibodytimecourseplusexp[1], a, b)* mersprobnoreinfectiontimecourse[day]
    
    mersprobreinfection.append(tempValue)      
    
import matplotlib.pyplot as plt7
plt7.plot(mersprobreinfection)
plt7.ylim(0, 0.0003)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('MERS') 

import matplotlib.pyplot as plt7

f = plt7.figure()

plt7.plot(mersprobreinfection)
plt7.ylim(0, 0.0003)
plt7.xlabel("Day")
plt7.ylabel("Probability of reinfection")
plt7.title('MERS') 
f.savefig("MERS_choe_PrInf-by-Day_" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-2

##### Input Data: Average peak-normalized ELIZA ODs for SARS-CoV-2 S IgG antibodies from Gudbjartsson et al. 2020, 
##### "Humoral Immune Response to SARS-CoV-2 in Iceland", New England Journal of Medicine.

###### days : 34, 48, 70, 94, 109
###### OD : 1.6, 1.55, 1.43, 1.19, 1.22

sarscov2data = [[1, 34, "SARSCoV2", 1.6/1.6, True, "NA", "NA"],
                [2, 48, "SARSCoV2", 1.55/1.6, False, (48 - 34), ((1.55 - 1.6)/1.6) /(48 - 34)],
                [3, 70, "SARSCoV2", 1.43/1.6, False, (70 - 48), ((1.43 - 1.55)/1.6) /(70 - 48)], 
                [4, 94, "SARSCoV2", 1.19/1.6, False, (94 - 70), ((1.19 - 1.43)/1.6)/(94 - 70)]]         
                
## SARS-CoV-2 Waning of Antibody OD

sarscov2datalength = len(sarscov2data)

for i in range(0,len(sarscov2data)):
    for j in range(0,len(sarscov2data[i])):
        print(sarscov2data[i][j])
        
aod = 0.75
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
    aod = aod + 0.05

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

while ( day < len(sarscov2paddedmeanwaning) and sarscov2antibodytimecourse[day] >= 0.75):
    print(day)
    print(sarscov2antibodytimecourse[day])
    print(sarscov2antibodytimecourse[day])
    print("sarscov2paddedmeanwaning")
    print(sarscov2paddedmeanwaning)
    print("sarscov2antibodytimecourse[day]")
    print(sarscov2antibodytimecourse[day])
    print("inerpolate of sarscov2paddedmeanwaning")
    #print(y_interpolation(sarscov2antibodytimecourse[day]))
    print("Added value")
    #test = sarscov2antibodytimecourse[day] + y_interpolation(sarscov2antibodytimecourse[day])
    #Not working
    #sarscov2antibodytimecourse.append(sarscov2antibodytimecourse[day] + y_interpolation(sarscov2antibodytimecourse[day]))
    sarscov2antibodytimecourse.append(sarscov2antibodytimecourse[day] + 0)
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
df4.to_csv('SARS-CoV-2-Antibody-Time-Course' + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".csv", index=False)

sarscov2baseline = 0.1298892

###This baseline peak-normalized N IgG antibody level for SARS-CoV-2 comes from the ancestral and descendent states analysis that used the baselines for the human "seasonal" coronaviruses to estimate the baselines for the zoonotic coronaviruses.

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

sarscov2halflife = np.log(2) / 0.0042668219735056195

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
maxDays = 7
while day < maxDays:
    day = day + 1
    print("day")
    print(day)
    exponentValueBaseLine = -0.0042668219735056195
    #Not working
    #tempValue = sarscov2baseline + (sarscov2antibodytimecourseplusexp[day-1] - sarscov2baseline) * math.exp(exponentValueBaseLine)
    tempValue = sarscov2baseline + (sarscov2antibodytimecourseplusexp[day] - sarscov2baseline) * math.exp(exponentValueBaseLine)    
    print('tempValue')
    print(tempValue)
    sarscov2antibodytimecourseplusexp.append(tempValue)
    
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

f.savefig("SARS-CoV-2_nOD-by-Day" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

##### SARS-CoV-2 Probability of Infection | Antibody OD

import math

def sarscov2probinfgivenaod(aod):
    exponentValue = 5.148993 + (10.89892 * aod)
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
plt11.ylim(0,0.006)
plt11.xlim(0.0, 1.0)
plt11.xlabel("Peak - normalized IgG antibody OD")
plt11.ylabel("Daily probability of infection")
plt11.title('SARS-CoV-2')
plt11.show()

import matplotlib.pyplot as plt12

f12 = plt12.figure()
plt12.plot(sarscov2probinfgivenaodList)

plt12.ylim(0,0.006)
plt12.xlim(0.0, 1.0)
plt12.xlabel("Peak - normalized IgG antibody OD")
plt12.ylabel("Daily probability of infection")
plt12.title('SARS-CoV-2')
plt12.show()

f12.savefig("SARS-CoV-2_PrInf-by-nOD" + str(datetime.now().strftime('%Y_%m_%d_%H_%M_%S')) + ".pdf", bbox_inches='tight')

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

from matplotlib import pyplot as plt13

sarscov2probinfList1 = populatesarscov1probinfList(-5.148993, -10.89892)

plt13.plot(sarscov2probinfList1, color='red')
plt13.ylim(0,0.008)
plt13.xlim(0.0, 1.0)
plt13.xlabel("Peak-normalized IgG antibody OD")

sarscov2probinfList2 = populatesarscov1probinfList(-4.229918, -15.59699)
sarscov2probinfList3 = populatesarscov1probinfList(-5.659934, -7.765139)
sarscov2probinfList4 = populatesarscov1probinfList(-5.007058, -11.626536)
sarscov2probinfList5 = populatesarscov1probinfList(-5.040719, -11.20207)

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
while day < 5:
    
    if(day < len(sarscov2antibodytimecourse)):
        sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day + 1])))
    else:
        #sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1392])))
        sarscov2probnoreinfectiontimecourse.append(sarscov2probnoreinfectiontimecourse[day] * (1 - sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[7])))
            
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
while day < 6:
    day = day + 1  
        
    #if day < 1392:
    if day < 6:
        tempValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[day]) * sarscov2probnoreinfectiontimecourse[day]
    else:       
        #empValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[1392])  * sarscov2probnoreinfectiontimecourse[day]
        tempValue = sarscov2probinfgivenaod(sarscov2antibodytimecourseplusexp[0])  * sarscov2probnoreinfectiontimecourse[day]
    
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

