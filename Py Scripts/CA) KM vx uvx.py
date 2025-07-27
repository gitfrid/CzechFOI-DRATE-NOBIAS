import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
import plotly.graph_objects as go
import os  # Used for extracting input filename


# Kaplan-Meier Survival Analysis: Vaccinated vs Unvaccinated

# This script performs a survival analysis using the Kaplan-Meier estimator
# on a dataset of individuals with vaccination and death dates. It compares
# survival between vaccinated and unvaccinated individuals, optionally filtered
# by age, and exports an interactive Plotly HTML graph.

#  - Input: CSV with birth year, death date, and up to 7 dose dates per person.
#  - Output: HTML file with Kaplan-Meier survival curves (total, vaccinated, unvaccinated).
#  - Requirements: pandas, numpy, lifelines, plotly, os


# === Constants ===

# Choose one of the input/output CSV file pairs by uncommenting as needed:

#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\CA) KM vx uvx\CA-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70.html"

INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\CA) KM vx uvx\CA) real data Vesely_106_202403141131_AG70.html"




# Configuration values
START_DATE = pd.Timestamp('2020-01-01')  # Reference date for day number conversion
MAX_AGE = 113                            # Max valid age
REFERENCE_YEAR = 2023                   # Used to calculate age from birth year

# === Age Filter ===
AGE_SELECTED = [70]  # Filter specific ages; use [] to include all ages

# === Load and Prepare Data ===
# Define dose date columns
dose_date_cols = [f'Datum_{i}' for i in range(1, 8)]
needed_cols = ['Rok_narozeni', 'DatumUmrti'] + dose_date_cols

# Load data with selected columns and parse date columns
df = pd.read_csv(
    INPUT_CSV,
    usecols=needed_cols,
    parse_dates=['DatumUmrti'] + dose_date_cols,
    dayfirst=False,
    low_memory=False
)

# Normalize column names
df.columns = [col.strip().lower() for col in df.columns]
dose_date_cols_lower = [col.lower() for col in dose_date_cols]

# Compute age from birth year
df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']

# Filter valid ages
df = df[df['age'].between(0, MAX_AGE)].copy()

# Apply age filter if any specific age(s) are selected
if AGE_SELECTED:
    df = df[df['age'].isin(AGE_SELECTED)]

# Function to convert dates to day numbers since START_DATE
def to_day_number(date_series):
    return (date_series - START_DATE).dt.days

# Convert death date and dose dates to day numbers
df['death_day'] = to_day_number(df['datumumrti'])
for col in dose_date_cols_lower:
    df[col + '_day'] = to_day_number(df[col])

# Compute the earliest dose day per person
df['first_dose_day'] = df[[col + '_day' for col in dose_date_cols_lower]].min(axis=1, skipna=True)

# Determine if person received any dose
df['has_any_dose'] = df[[col + '_day' for col in dose_date_cols_lower]].notna().any(axis=1)

# === Add censoring information ===
df['censor_day'] = df['death_day'].isna()

# Assign maximum observed day + 1 to censored observations
df['death_day'] = df['death_day'].fillna(df['death_day'].max() + 1)

# Event indicator: 1 = death observed, 0 = censored
df['event'] = (~df['censor_day']).astype(int)

# === Group assignment: vaccinated vs unvaccinated ===
df['group'] = 'uvx'  # Default: unvaccinated
df.loc[df['has_any_dose'], 'group'] = 'vx'  # Mark as vaccinated if any dose
df['group'] = df['group'].astype('category')  # Optimize memory usage

# === Debug prints ===
print("Data shape after filtering by age:", df.shape)
print("Group counts:\n", df['group'].value_counts())
print("Death day stats:")
print(df['death_day'].describe())
print("Event counts:")
print(df['event'].value_counts())
print("Any NaNs in death_day?", df['death_day'].isna().sum())

# === Fit Kaplan-Meier curves for each group ===
kmf_total = KaplanMeierFitter()
kmf_vx = KaplanMeierFitter()
kmf_uvx = KaplanMeierFitter()

# Extract survival durations and event indicators
T_total = df['death_day']
E_total = df['event']

T_vx = df[df['group'] == 'vx']['death_day']
E_vx = df[df['group'] == 'vx']['event']

T_uvx = df[df['group'] == 'uvx']['death_day']
E_uvx = df[df['group'] == 'uvx']['event']

# Fit KM models
kmf_total.fit(T_total, event_observed=E_total, label='Total')
kmf_vx.fit(T_vx, event_observed=E_vx, label='Vaccinated')
kmf_uvx.fit(T_uvx, event_observed=E_uvx, label='Unvaccinated')

# Get base name of input file for plot subtitle
input_filename = os.path.basename(INPUT_CSV)

# === Create Plotly figure ===
fig = go.Figure()

# Add KM survival curves to figure
for kmf in [kmf_total, kmf_vx, kmf_uvx]:
    fig.add_trace(go.Scatter(
        x=kmf.survival_function_.index,
        y=kmf.survival_function_[kmf._label],
        mode='lines',
        name=kmf._label
    ))

# Update layout with titles and labels
fig.update_layout(
    title=f'Kaplan-Meier Survival Curves: Total vs Vaccinated vs Unvaccinated AGE:{AGE_SELECTED}<br><sub>Input CSV: {input_filename}</sub>',
    xaxis_title='Days Since Jan 1, 2020',
    yaxis_title='Survival Probability',
    template='plotly_white'
)

# Export to HTML file
fig.write_html(OUTPUT_HTML)
print(f"Plot saved to: {OUTPUT_HTML}")
