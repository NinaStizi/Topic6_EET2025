"""
 Topic 6 - Effect of environmental regulation on Water Pollution and Infant Mortality in India

 --
 Extensions of Hanna and Greenstone, American Economic Review, 2014

    1. DiD Extension
    2. DDD Design
    3. Robustness Checks

"""
#***************************************#
#********** GENERAL SET UP *************#
#***************************************#

# Packages
import os  
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf


# Automating the setup of the right directory
script_dir = os.path.dirname(os.path.abspath(__file__))
mypath = os.path.join(os.path.dirname(script_dir))
os.chdir(mypath)
print(os.getcwd())

#***************************************#
#********** DATA CLEANING **************#
#***************************************#

# Data Loading and Cleaning
cities_df = pd.read_stata("data/india_waters_cityyear.dta")

cities_df = cities_df[["year", "river", "state", "city", "district", "nrcp", "bod", "fcoli", "lnfcoli", "do", "corruption_score05", "corruption_level07", "top5_corrupt", "bottom5_corrupt", "pop_urban", "lit_urban", 'povgap', 'total_industries']]

# First clean : 36 cities out of 425 are NAs, so we drop them
cities_df = cities_df.dropna(how = 'all')
cities_df = cities_df.reset_index(drop = True)

# Since we will be working on corruption data, we drop the 11 cities out of 425 for which there are no corruption level / corruption scores
cities_df = cities_df[~cities_df.corruption_level07.isin([''])]


#***************************************#
#*********** 1. DiD EXTENSION **********#
#***************************************#

#***************************************#
#*********** 2. DDD DESIGN *************#
#***************************************#

# Creation of variables of interest for regression:

cities_df['corruption'] = np.where(
    (cities_df['corruption_level07'] == 'Moderate') & (cities_df['corruption_score05'] <= 480),
    0,
    1
)

cities_df["post"] = (cities_df["year"] >= 2000).astype(int)

cities_df["city_fe"] = cities_df["city"]
cities_df["year_fe"] = cities_df["year"].astype(str)  # 

cities_df["nrcp_post"] = cities_df["nrcp"] * cities_df["post"]
cities_df["nrcp_post_corruption"] = cities_df["nrcp_post"] * cities_df["corruption"]

model_df = cities_df[[
    'bod', 'do','fcoli', 'nrcp', 'post', 'corruption', 
    'pop_urban', 'lit_urban', 'povgap', 'total_industries', 
    'city'
]].dropna()

model_df['nrcp_post'] = model_df['nrcp'] * model_df['post']
model_df['nrcp_post_corruption'] = model_df['nrcp_post'] * model_df['corruption']


# BOD
model_bod_ddd = smf.ols(
    formula="""
        bod ~ nrcp + post + corruption + nrcp:post + 
        nrcp:corruption + post:corruption + 
        nrcp_post_corruption + 
        pop_urban + lit_urban + povgap + total_industries
    """,
    data=model_df
).fit(cov_type='cluster', cov_kwds={'groups': model_df['city']})

# DO
model_do_ddd = smf.ols(
    formula="""
        do ~ nrcp + post + corruption + nrcp:post + 
        nrcp:corruption + post:corruption + 
        nrcp_post_corruption + 
        pop_urban + lit_urban + povgap + total_industries
    """,
    data=model_df
).fit(cov_type='cluster', cov_kwds={'groups': model_df['city']})

# FColi
model_fcoli_ddd = smf.ols(
    formula="""
        fcoli ~ nrcp + post + corruption + nrcp:post + 
        nrcp:corruption + post:corruption + 
        nrcp_post_corruption + 
        pop_urban + lit_urban + povgap + total_industries
    """,
    data=model_df
).fit(cov_type='cluster', cov_kwds={'groups': model_df['city']})

# Printing the summaries
print(model_bod_ddd.summary())
print(model_do_ddd.summary())
print(model_fcoli_ddd.summary())

#***************************************#
#********* 3. ROBUSTNESS CHECKS ********#
#***************************************#