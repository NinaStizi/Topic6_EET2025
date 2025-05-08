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
import pandas as pd
import numpy as np
import os  

# Automating the setup of the right directory
script_dir = os.path.dirname(os.path.abspath(__file__))
mypath = os.path.join(os.path.dirname(script_dir))
os.chdir(mypath)
print(os.getcwd())

#**********************************************************************************************#
#**** The role of corruption and other city characteristics on the failure of NRCP policy *****#
#**********************************************************************************************#

# Data Loading
cities_df = pd.read_stata("data/india_waters_cityyear.dta")
