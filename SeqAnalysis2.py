import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

#Open pickled output from SeqCollection.py
grids = pickle.load(open("grids.p",'rb'))

#Want to convert pandas datatables into numpy arrays so they work with legacy code for heatmaps
def norm(matrix):
    matrix = (matrix.values)
    norm = np.sum(matrix)
    print(norm)
    matrix = matrix/float(norm)
    return(np.log(matrix))

#Rename and reformat grids from pickle file for visulization
ts1 = norm(grids['CobaltT0R1'])
ts2 = norm(grids['CobaltT1R1'])
ts3 = norm(grids['CobaltT2R1'])
ts4 = norm(grids['CobaltT0R2'])
ts5 = norm(grids['CobaltT1R2'])
ts6 = norm(grids['CobaltT2R2'])
ts7 = norm(grids['ControlT0R1'])
ts8 = norm(grids['ControlT1R1'])
ts9 = norm(grids['ControlT2R1'])

#Initialize slope-containing grids
slope_grid1 = np.zeros((21,76))
slope_grid2 = np.zeros((21,76))

#loop over grid
for i in range(21):
    for j in range(76):
        x1 = [0.0, 2.129283017, 4.794863978] #Treated population doublings at each time point
        x2 = [0.0, 2.022367813, 3.50041511] #Control population doublings
        y1 = [ts4[i][j], ts5[i][j], ts6[i][j]] #log frequency of amino acid i, position j at t0, t1, and t2 for treated cells
        y2 = [ts7[i][j], ts8[i][j], ts9[i][j]] #log frequency of amino acid i, position j at t0, t1, and t2 for control cells

        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1,y1) # Calculate change in log frequency / generation in treated
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2,y2) # Same for control cells
        
        slope_grid1[i][j] = slope1 # assign slope to appropriate grid entry
        slope_grid2[i][j] = slope2 # assign slope to appropriate grid entry


#Define amino acids in order for convenient plotting later on:
labels = ['STOP', 'W', 'F', 'Y', 'L', 'I', 'M', 'V', 'C', 'A', 'G', 'P', 'S',
'T', 'N', 'Q', 'H', 'R', 'K', 'D', 'E']

##############################################################
# Figure 1: Compare slopes for individual cells in the grid for control and treated cells
##############################################################
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

##############################################################
# Figure 2: Plot control slopes v treated slopes to look for correlation and outlying points
##############################################################
plt.figure(figsize = (6,6), facecolor = 'white')
x = np.linspace(-1,1,10)
plt.plot(x,x,ls = '--', color = 'black', lw = 2, alpha = 0.5) # y=x reference line
plt.scatter(slope_grid1,slope_grid2, marker = 'o', color = 'blue', alpha = 0.3)
plt.axis([-1,1,-1,1])
plt.xlabel('Control Growth Rate')
plt.ylabel('Treated Growth Rate')
plt.grid()
plt.show()
