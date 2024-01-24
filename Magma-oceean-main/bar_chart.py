import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd


Features=['T(K)', 'Na2O', 'MgO', 'Al2O3', 'SiO2', 'K2O', 'CaO', 'TiO2', 'MnO', 'FeO', 'Fe2O3']
Feature_importance = [0.18663728, 0.00853247, 0.06739053, 0.02869617, 0.01628545,
       0.02380994, 0.02506919, 0.00564358, 0.00450428, 0.21438232,
       0.4190488 ]

plt.bar(range(len(Feature_importance)),Feature_importance, width = 1, tick_label = Features, edgecolor= 'black',color='#eb4634',alpha=0.8)
plt.xlabel("Features")
plt.ylabel("Feature Importance")
plt.title("Feature Importance from Random Forest Regression")
plt.show()
