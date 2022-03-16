# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:06:28 2022

@author: rossc
"""

import matplotlib.pyplot as plt
import pandas as pd

   
dataframe = pd.read_csv("Training.csv")
print(dataframe)

x = dataframe.Epoch
y = dataframe.EnergyRMSE

plt.plot(x, y, color="b", lw=1.5)
plt.xlabel("Epoch")
plt.ylabel("Energy RMSE (kcal/mol)")
plt.title("N = 500", fontsize=20 )
plt.xticks()
plt.yticks()
plt.savefig("EnergyRMSEvEpoch.png", dpi=600)
plt.show()