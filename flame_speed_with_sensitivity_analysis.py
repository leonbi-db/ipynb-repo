# Databricks notebook source
# MAGIC %md
# MAGIC # Flame Speed with Sensitivity Analysis

# COMMAND ----------

# MAGIC %md
# MAGIC In this example we simulate a freely-propagating, adiabatic, 1-D flame and
# MAGIC * Calculate its laminar burning velocity
# MAGIC * Perform a sensitivity analysis of its kinetics
# MAGIC 
# MAGIC The figure below illustrates the setup, in a flame-fixed co-ordinate system. The reactants enter with density $\rho_{u}$, temperature $T_{u}$ and speed $S_{u}$. The products exit the flame at speed $S_{b}$, density $\rho_{b}$ and temperature $T_{b}$.

# COMMAND ----------

# MAGIC %md
# MAGIC <img src="images/flameSpeed.png" alt="Freely Propagating Flame" style="width: 300px;"/>

# COMMAND ----------

# MAGIC %md
# MAGIC ### Import Modules

# COMMAND ----------

from __future__ import print_function
from __future__ import division

import cantera as ct
import numpy as np

print("Running Cantera Version: " + str(ct.__version__))

# COMMAND ----------

# Import plotting modules and define plotting preference
%matplotlib notebook
import matplotlib.pylab as plt

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (8,6)

# Get the best of both ggplot and seaborn
plt.style.use('ggplot')
plt.style.use('seaborn-deep')

plt.rcParams['figure.autolayout'] = True

# Import Pandas for DataFrames
import pandas as pd

# COMMAND ----------

# MAGIC %md
# MAGIC ### Define the reactant conditions, gas mixture and kinetic mechanism associated with the gas

# COMMAND ----------

#Inlet Temperature in Kelvin and Inlet Pressure in Pascals
#In this case we are setting the inlet T and P to room temperature conditions
To = 300
Po = 101325

#Define the gas-mixutre and kinetics
#In this case, we are choosing a GRI3.0 gas
gas = ct.Solution('gri30.cti')

# Create a stoichiometric CH4/Air premixed mixture 
gas.set_equivalence_ratio(1.0, 'CH4', {'O2':1.0, 'N2':3.76})
gas.TP = To, Po

# COMMAND ----------

# MAGIC %md
# MAGIC ### Define flame simulation conditions

# COMMAND ----------

# Domain width in metres
width = 0.014

# Create the flame object
flame = ct.FreeFlame(gas, width=width)

# Define tolerances for the solver
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

# Define logging level
loglevel = 1

# COMMAND ----------

# MAGIC %md
# MAGIC ### Solve

# COMMAND ----------

flame.solve(loglevel=loglevel, auto=True)
Su0 = flame.u[0]
print("Flame Speed is: {:.2f} cm/s".format(Su0*100))

# Note that the variable Su0 will also be used downsteam in the sensitivity analysis

# COMMAND ----------

# MAGIC %md
# MAGIC ### Plot figures
# MAGIC 
# MAGIC Check and see if all has gone well. Plot temperature and species fractions to see
# MAGIC 
# MAGIC #### Temperature Plot

# COMMAND ----------

plt.figure()

plt.plot(flame.grid*100, flame.T, '-o')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (K)');

# COMMAND ----------

# MAGIC %md
# MAGIC #### Major species' plot

# COMMAND ----------

# MAGIC %md
# MAGIC To plot species, we first have to identify the index of the species in the array
# MAGIC For this, cut & paste the following lines and run in a new cell to get the index
# MAGIC 
# MAGIC     for i, specie in enumerate(gas.species()):
# MAGIC         print(str(i) + '. ' + str(specie))

# COMMAND ----------

# Extract concentration data
X_CH4 = flame.X[13]
X_CO2 = flame.X[15]
X_H2O = flame.X[5]

plt.figure()

plt.plot(flame.grid*100, X_CH4, '-o', label=r'$CH_{4}$')
plt.plot(flame.grid*100, X_CO2, '-s', label=r'$CO_{2}$')
plt.plot(flame.grid*100, X_H2O, '-<', label=r'$H_{2}O$')

plt.legend(loc=2)
plt.xlabel('Distance (cm)')
plt.ylabel('MoleFractions');

# COMMAND ----------

# MAGIC %md
# MAGIC ## Sensitivity Analysis
# MAGIC 
# MAGIC See which reactions effect the flame speed the most

# COMMAND ----------

# Create a dataframe to store sensitivity-analysis data
sensitivities = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))

# COMMAND ----------

# MAGIC %md
# MAGIC ### Compute sensitivities

# COMMAND ----------

# Set the value of the perturbation
dk = 1e-2

# Create an empty column to store the sensitivities data
sensitivities["baseCase"] = ""

# COMMAND ----------

for m in range(gas.n_reactions):
    gas.set_multiplier(1.0) # reset all multipliers                                                                     
    gas.set_multiplier(1+dk, m) # perturb reaction m   
    
    # Always force loglevel=0 for this
    # Make sure the grid is not refined, otherwise it won't strictly 
    # be a small perturbation analysis
    flame.solve(loglevel=0, refine_grid=False)
    
    # The new flame speed
    Su = flame.u[0]
    
    sensitivities["baseCase"][m] = (Su-Su0)/(Su0*dk)

# This step is essential, otherwise the mechanism will have been altered
gas.set_multiplier(1.0)

# COMMAND ----------

sensitivities.head()

# COMMAND ----------

# MAGIC %md
# MAGIC ### Make plots

# COMMAND ----------

# Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
# to see only the top few
threshold = 0.03

firstColumn = sensitivities.columns[0]

# For plotting, collect only those steps that are above the threshold
# Otherwise, the y-axis gets crowded and illegible
sensitivitiesSubset = sensitivities[sensitivities[firstColumn].abs() > threshold]
indicesMeetingThreshold = sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for GRI 3.0",
                                                          legend=None)
plt.gca().invert_yaxis()

plt.rcParams.update({'axes.labelsize': 20})
plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');

# Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
# plt.savefig('sensitivityPlot', dpi=300)
