"""
 Topic 6 - Effect of environmental regulation on Water Pollution and Infant Mortality in India

 --
 Extensions of Hanna and Greenstone, American Economic Review, 2014

    1. DiD Extension: from the model of the authors, we add a specification that considers the level of corruption in the pollution level
    2. Event Study Analysis, using Sun and Abraham (2021) method

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
import statsmodels.api as sm


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

cities_df = cities_df[["year", "river", "state", "city", "district", "nrcp", "bod", "lnfcoli", "do", "corruption_level07", "corruption_score05", "pop_urban", "lit_urban", 'povgap', 'pce', 'total_industries']]

# First clean : 36 cities out of 425 are NAs, so we drop them
cities_df = cities_df.dropna(how = 'all')
# Some cities have no river... we drop them too
cities_df = cities_df[~cities_df.river.isin([''])]


# Since we will be working on corruption data, we drop the 11 cities out of 425 for which there are no corruption level / corruption scores
cities_df = cities_df[~cities_df.corruption_level07.isin([''])]

# There are 3 cities for which NRCP is not 1 or 0 we drop them
cities_df = cities_df[cities_df.nrcp.isin([0,1])]

# We order the dataset by state, river, city and year
cities_df = cities_df.sort_values(by = ["state", "river", "city", "year"], ascending = True)
cities_df = cities_df.reset_index(drop = True)

# We create the cohorts:
cities_df['nrcp_year'] = cities_df.apply(
    lambda row: row['year'] if row['nrcp'] == 1 else pd.NA, axis=1
)
cities_df['nrcp_implementation'] = (
    cities_df.groupby('city')['nrcp_year']
    .transform('min')
)

unique_implementations = (
    cities_df.dropna(subset=['nrcp_implementation'])
    .drop_duplicates(subset=['city'])[['city', 'nrcp_implementation']]
)

unique_implementations['cohort'] = unique_implementations['nrcp_implementation'].rank(method='dense').astype(int)

cities_df = cities_df.merge(unique_implementations[['city', 'cohort']], on='city', how='left')

cities_df['cohort'] = cities_df['cohort'].fillna('control')
cities_df['nrcp_implementation'] = cities_df['nrcp_implementation'].fillna('control')



# We create the list of variables that we will use :
pollutants = ["bod", "lnfcoli", "do"]

controls = ["pop_urban", "lit_urban", 'povgap', 'pce', 'total_industries']

corruption_levels = cities_df.corruption_level07.unique()


#***************************************#
#*********** 1. DiD EXTENSION **********#
#***************************************#

for pollutant in pollutants:
    formula = f"""
        {pollutant} ~ nrcp + corruption_level07 + nrcp * C(corruption_level07) + {' + '.join(controls)} + C(city) + C(year)
        
    """

    # Drop missing data used in the current regression
    model_data = cities_df[[pollutant, 'nrcp', 'city', 'year', 'corruption_level07'] + controls].dropna()

    # Reattach city for clustering
    model_data['city'] = model_data['city']

    model = smf.ols(
        formula=formula,
        data=model_data
    ).fit(cov_type='cluster', cov_kwds={'groups': model_data['city']})

    # Filter coefficients (keep only the most relevant)
    print(f"\n=== Results for {pollutant.upper()} ===")
    print(model.summary())

    # Save LaTeX summary
    latex_code = model.summary().as_latex()
    with open(f"results/corruption/{pollutant}_DiD_regression.tex", "w") as f:
        f.write(latex_code)

    print(f"Saved LaTeX regression table for {pollutant} ✅")


#***************************************#
#******** 2. EVENT STUDY ANALYSIS ******#
#***************************************#

# We introduce the treatment year per city
cities_df['treatment_year'] = cities_df.groupby('city')['year'].transform(
    lambda x: cities_df.loc[x.index, 'nrcp'].eq(1).replace(False, pd.NA) * cities_df.loc[x.index, 'year']
)
cities_df['first_treatment'] = cities_df.groupby('city')['treatment_year'].transform('min')

# We create event time (relative year to treatment)
cities_df['event_time'] = cities_df['year'] - cities_df['first_treatment']

# We create safe event-time labels again
event_time_labels = {}
for t in range(-5, 6):
    if t == -1: continue  # base year
    label = f"event_m{abs(t)}" if t < 0 else f"event_p{t}"
    cities_df[label] = (cities_df['event_time'] == t).astype(int)
    event_time_labels[t] = label

# We create interaction terms: corruption * event dummies
for label in event_time_labels.values():
    cities_df[f"{label}_x_corr"] = cities_df[label] * cities_df["corruption_score05"]

# We run the Sun and Abraham style regression for each pollutant
for pollutant in pollutants:
    print(f"\n=== Corruption Effect on {pollutant.upper()} Around NRCP ===")

    event_corr_terms = [f"{label}_x_corr" for label in event_time_labels.values()]
    formula = f"""
        {pollutant} ~ {' + '.join(event_corr_terms + controls)} + C(city) + C(year)
    """

    model_data = cities_df[[pollutant, 'city', 'year', 'corruption_score05'] + event_corr_terms + controls].dropna()

    model = smf.ols(formula=formula, data=model_data).fit(cov_type="cluster", cov_kwds={"groups": model_data["city"]})

    # Tables to show the marginal effect of corruption by event time
    coef_table = model.summary2().tables[1].loc[event_corr_terms]
    print(coef_table)

    table_latex = coef_table.to_latex()
    with open(f"results/corruption/{pollutant}_corruption_effect_SunAbraham.tex", "w") as f:
        f.write(table_latex)

    # Plots to shows marginal effect of corruption on pollution over time
    plt.figure(figsize=(8, 5))
    years = list(event_time_labels.keys())
    coef = model.params[event_corr_terms]
    err = model.bse[event_corr_terms]

    plt.errorbar(years, coef, yerr=1.96 * err, fmt='o-', capsize=4, color='darkgreen')
    plt.axhline(0, linestyle='--', color='black')
    plt.axvline(0, linestyle=':', color='red', label='NRCP Implementation')
    plt.title(f"Effect of Corruption on {pollutant.upper()} by Event Year")
    plt.xlabel("Years since NRCP treatment")
    plt.ylabel("Marginal Effect of Corruption on Pollutant")
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"results/corruption/{pollutant}_corruption_effect_SunAbraham.pdf")
    plt.show()
    plt.close()