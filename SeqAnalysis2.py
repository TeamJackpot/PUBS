import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

grids = pickle.load(open("grids.p",'rb'))

def norm(matrix):
    matrix = (matrix.values)
    norm = np.sum(matrix)
    print(norm)
    matrix = matrix/float(norm)
    return(np.log(matrix))

ts1 = norm(grids['CobaltT0R1'])
ts2 = norm(grids['CobaltT1R1'])
ts3 = norm(grids['CobaltT2R1'])
ts4 = norm(grids['CobaltT0R2'])
ts5 = norm(grids['CobaltT1R2'])
ts6 = norm(grids['CobaltT2R2'])
ts7 = norm(grids['ControlT0R1'])
ts8 = norm(grids['ControlT1R1'])
ts9 = norm(grids['ControlT2R1'])

slope_grid1 = np.zeros((21,76))
slope_grid2 = np.zeros((21,76))

for i in range(21):
    for j in range(76):
        x1 = [0.0, 2.129283017, 4.794863978]
        x2 = [0.0, 2.022367813, 3.50041511]
        y1 = [ts4[i][j], ts5[i][j], ts6[i][j]]
        y2 = [ts7[i][j], ts8[i][j], ts9[i][j]]
        #y1 = [ts4[i][j], ts5[i][j], ts6[i][j]]
        y2 = [ts7[i][j], ts8[i][j], ts9[i][j]]
        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1,y1)
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2,y2)
        slope_grid1[i][j] = slope1
        slope_grid2[i][j] = slope2

plt.figure(figsize = (6,6), facecolor = 'white')
plt.scatter(slope_grid1,slope_grid2, marker = 'o', color = 'blue', alpha = 0.3)
x = np.linspace(-1,1,10)
plt.plot(x,x,ls = '--', color = 'black', lw = 2, alpha = 0.5)
plt.axis([-1,1,-1,1])
plt.xlabel('Control Growth Rate')
plt.ylabel('Treated Growth Rate')
plt.grid()
plt.show()

labels = ['STOP', 'W', 'F', 'Y', 'L', 'I', 'M', 'V', 'C', 'A', 'G', 'P', 'S',
'T', 'N', 'Q', 'H', 'R', 'K', 'D', 'E']

plt.figure(figsize = (6,6), facecolor = 'white')

plt.subplot(211)
plt.imshow(slope_grid1, cmap = 'bwr', interpolation = 'nearest',aspect = 1, vmin=-1,vmax=1)
plt.xlabel("Position")
plt.ylabel("Amino Acid (Coded)")
plt.title("Control Sample Slope")
plt.yticks(range(21),labels, rotation = 'horizontal')
plt.colorbar(label = "Change in (Log(Frequency) over Time)")

plt.subplot(212)
plt.imshow(slope_grid2, cmap = 'bwr', interpolation = 'nearest',aspect = 1, vmin=-1,vmax=1)
plt.xlabel("Position")
plt.ylabel("Amino Acid (Coded)")
plt.title("Treated Sample Slope")
plt.yticks(range(21),labels, rotation = 'horizontal')
plt.colorbar(label = "Change in (Log(Frequency) over Time)")
plt.show()

plt.figure(figsize = (6,6), facecolor = 'white')
plt.imshow(np.log(ts3), cmap = 'Blues', interpolation = 'nearest',aspect = 1)
plt.xlabel("Position")
plt.ylabel("Amino Acid (Coded)")
plt.title("Control Counts at T = 1")
plt.yticks(range(21),labels, rotation = 'horizontal')
plt.colorbar(label = "Log Counts")

'''
plt.subplot(324)
plt.imshow(ts2, cmap = 'Blues', interpolation = 'nearest',aspect = 1)
plt.xlabel("Position")
plt.ylabel("Amino Acid (Coded)")
plt.title("Treated Counts at T = 1")
plt.yticks(range(21),labels, rotation = 'horizontal')
#plt.colorbar(label = "Fitness Index")

plt.subplot(325)
plt.imshow(ts9, cmap = 'Blues', interpolation = 'nearest',aspect = 1)
plt.xlabel("Position")
plt.ylabel("Amino Acid (Coded)")
plt.title("Control Counts at T = 2")
plt.yticks(range(21),labels, rotation = 'horizontal')
#plt.colorbar(label = "Fitness Index")

plt.subplot(326)
plt.imshow(ts3, cmap = 'Blues', interpolation = 'nearest',aspect = 1)
plt.xlabel("Position")
plt.ylabel("Amino Acid (Coded)")
plt.title("Treated Counts at T = 2")
plt.yticks(range(21),labels, rotation = 'horizontal')
#plt.colorbar(label = "Fitness Index")
'''
