"""
 Topic 6 - Effect of environmental regulation on Water Pollution and Infant Mortality in India

 --
 Extensions of Hanna and Greenstone, American Economic Review, 2014

    1. The role of corruption and other city characteristics

"""
#**********************************************************************************************#
#**************** General SetUp : Environment, Working Directory, Packages ********************#
#**********************************************************************************************#

# Packages
import os  
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Automating the setup of the right directory
script_dir = os.path.dirname(os.path.abspath(__file__))
mypath = os.path.join(os.path.dirname(script_dir))
os.chdir(mypath)
print(os.getcwd())

#**********************************************************************************************#
#**** The role of corruption and other city characteristics on the failure of NRCP policy *****#
#**********************************************************************************************#

# Data Loading and Cleaning
cities_df = pd.read_stata("data/india_waters_cityyear.dta")

cities_df = cities_df[["year", "river", "state", "city", "district", "nrcp", "bod", "fcoli", "lnfcoli", "do", "corruption_references", "pollution_references", "water_references", "corruption_score05", "corruption_level07", "top5_corrupt", "bottom5_corrupt" ]]

# Dummy: Visual observation
cities_per_corruptionlevel = cities_df.groupby(['year', 'corruption_level07'])[['bod', 'fcoli', 'lnfcoli', 'do']].mean()

corruption_levels_to_drop = [''] #droppting the undefined
cities_per_corruptionlevel = cities_per_corruptionlevel[~cities_per_corruptionlevel.index.get_level_values('corruption_level07').isin(corruption_levels_to_drop)]

## For BOD
plt.figure(figsize=(10, 6))

for corruption_level in cities_per_corruptionlevel.index.get_level_values('corruption_level07').unique():
    filtered_data = cities_per_corruptionlevel.xs(corruption_level, level='corruption_level07')
    plt.plot(filtered_data.index, filtered_data['bod'], label=f'{corruption_level}')

plt.title('Evolution of BOD Over Time by Corruption Level')
plt.xlabel('Year')
plt.ylabel('Mean BOD')
plt.legend(title='Corruption Level')
plt.grid(True)

plt.show()

## For FColi
plt.figure(figsize=(10, 6))

for corruption_level in cities_per_corruptionlevel.index.get_level_values('corruption_level07').unique():
    filtered_data = cities_per_corruptionlevel.xs(corruption_level, level='corruption_level07')
    plt.plot(filtered_data.index, filtered_data['fcoli'], label=f'{corruption_level}')

plt.title('Evolution of BOD Over Time by Corruption Level')
plt.xlabel('Year')
plt.ylabel('Mean BOD')
plt.legend(title='Corruption Level')
plt.grid(True)

plt.show()


