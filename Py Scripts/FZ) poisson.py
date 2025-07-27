import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import KaplanMeierFitter
import plotly.graph_objects as go
import sys

# === Constants and input ===
"""
Poisson Regression and Survival Analysis of Simulated or Real Death and Vaccination Data by Age

This script performs a Poisson regression and Kaplan-Meier survival analysis on individual-level 
death and vaccination data. It estimates incidence rate ratios (IRRs) for vaccinated vs. unvaccinated 
individuals and plots survival curves.

Key Steps:
1. Load and preprocess death and dose date data.
2. Expand each individual's timeline into daily person-day records.
3. Label each day with vaccination status and death occurrence.
4. Fit a Poisson regression model to estimate IRRs.
5. Plot Kaplan-Meier survival curves for vaccinated and unvaccinated groups.

Required:
- Input CSV with 'Rok_narozeni', 'DatumUmrti', and 'Datum_1' to 'Datum_7'
"""


# === Constants and input file paths ===
INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.TXT"

#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ) real data Vesely_106_202403141131_AG70 poisson.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ) real data Vesely_106_202403141131_AG70 poisson.TXT"

START_DATE = pd.Timestamp('2020-01-01')  # Reference date for day numbering
MAX_AGE = 113                            # Max age cutoff
REFERENCE_YEAR = 2023                    # Year used to calculate age from birth year

# === Logger that writes to both stdout and a log file ===
class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, data):
        for f in self.files:
            f.write(data)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()

# Redirect stdout/stderr to both console and log file
original_stdout = sys.stdout
original_stderr = sys.stderr
log_file_path = fr"{OUTPUT_TXT}"
log_file = open(log_file_path, "w", encoding="utf-8")
tee = Tee(sys.stdout, log_file)
sys.stdout = tee
sys.stderr = tee

# === Load and Prepare Data ===

# Specify columns needed: birth year, death date, and up to 7 dose dates
dose_date_cols = [f'Datum_{i}' for i in range(1, 8)]
needed_cols = ['Rok_narozeni', 'DatumUmrti'] + dose_date_cols

# Load selected columns with date parsing
df = pd.read_csv(
    INPUT_CSV,
    usecols=needed_cols,
    parse_dates=['DatumUmrti'] + dose_date_cols,
    dayfirst=False,
    low_memory=False
)

# Convert column names to lowercase for consistency
df.columns = [col.strip().lower() for col in df.columns]
dose_date_cols_lower = [col.lower() for col in dose_date_cols]

# Compute age from birth year
df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']

# Filter out invalid or extreme ages
df = df[df['age'].between(0, MAX_AGE)].copy()

# Helper function: convert datetime to day number since START_DATE
def to_day_number(date_series):
    return (date_series - START_DATE).dt.days

# Convert death and dose dates to day numbers
df['death_day'] = to_day_number(df['datumumrti'])
for col in dose_date_cols_lower:
    df[col + '_day'] = to_day_number(df[col])

# Find the earliest dose day and whether any dose was received
df['first_dose_day'] = df[[col + '_day' for col in dose_date_cols_lower]].min(axis=1, skipna=True)
df['has_any_dose'] = df[[col + '_day' for col in dose_date_cols_lower]].notna().any(axis=1)

# === Define Observation Period ===

# The end of measurement is defined by the latest death day
END_MEASURE = int(df['death_day'].dropna().max())
print(f"END_MEASURE (max death day): {END_MEASURE}")

# For censoring: alive = END_MEASURE, dead = death_day
df['end_day'] = df['death_day'].fillna(END_MEASURE)

# === Expand dataset to daily person-day records per individual ===

def expand_person_days(row):
    days_range = np.arange(0, int(row['end_day']) + 1)
    vacc_status = (days_range >= row['first_dose_day']) & row['has_any_dose']
    return pd.DataFrame({
        'age': row['age'],
        'day': days_range,
        'vaccinated': vacc_status.astype(int),
        'death': 0
    })

print("Expanding person-day data... This may take some time for large data.")

# Expand all individuals into person-day records
person_days_list = []
for idx, row in df.iterrows():
    person_days_list.append(expand_person_days(row))

person_days = pd.concat(person_days_list, ignore_index=True)

# Mark death occurrence on the corresponding person-day
death_rows = df.loc[df['death_day'].notna(), ['age', 'death_day']]
for _, death_row in death_rows.iterrows():
    mask = (person_days['age'] == death_row['age']) & (person_days['day'] == death_row['death_day'])
    person_days.loc[mask, 'death'] = 1

# === Aggregate Data for Poisson Regression ===

# Group by age, day, and vaccination status
agg = person_days.groupby(['age', 'day', 'vaccinated']).agg(
    deaths=('death', 'sum'),
    person_days=('day', 'count')
).reset_index()

# Poisson model requires offset = log(person-time), and we center age
agg['offset'] = np.log(agg['person_days'])
agg['age_c'] = agg['age'] - agg['age'].mean()

# Prepare design matrix for Poisson regression
X = sm.add_constant(agg[['vaccinated', 'age_c']])

# Fit Poisson regression: deaths ~ vaccinated + age_c + offset(log person-days)
model = sm.GLM(agg['deaths'], X, offset=agg['offset'], family=sm.families.Poisson())
result = model.fit()

# Output regression results
print(result.summary())

# Compute Incidence Rate Ratios (IRR) and 95% confidence intervals
params = result.params
conf = result.conf_int()
irr = np.exp(params)
irr_conf_lower = np.exp(conf[0])
irr_conf_upper = np.exp(conf[1])

print("\nIncidence Rate Ratios (IRRs):")
print(f"Intercept: {irr['const']:.3f}")
print(f"Vaccinated vs Unvaccinated: {irr['vaccinated']:.3f} (95% CI: {irr_conf_lower['vaccinated']:.3f} - {irr_conf_upper['vaccinated']:.3f})")

# === Kaplan-Meier Survival Plot ===

# Construct survival intervals for unvaccinated period
unvaccinated = df[['age', 'end_day', 'first_dose_day', 'death_day']].copy()
unvaccinated['start'] = 0
unvaccinated['stop'] = unvaccinated['first_dose_day'].fillna(unvaccinated['end_day'])
unvaccinated['event'] = (unvaccinated['death_day'] <= unvaccinated['stop']) & (unvaccinated['death_day'].notna())
unvaccinated['event'] = unvaccinated['event'].astype(int)
unvaccinated['group'] = 'Unvaccinated'
unvaccinated = unvaccinated[unvaccinated['stop'] > unvaccinated['start']]

# Construct survival intervals for vaccinated period
vaccinated = df[df['first_dose_day'].notna()][['age', 'end_day', 'first_dose_day', 'death_day']].copy()
vaccinated['start'] = vaccinated['first_dose_day']
vaccinated['stop'] = vaccinated['end_day']
vaccinated['event'] = (vaccinated['death_day'] >= vaccinated['start']) & (vaccinated['death_day'].notna())
vaccinated['event'] = vaccinated['event'].astype(int)
vaccinated['group'] = 'Vaccinated'
vaccinated = vaccinated[vaccinated['stop'] > vaccinated['start']]

# Combine and compute duration
km_data = pd.concat([unvaccinated, vaccinated], ignore_index=True)
km_data['duration'] = km_data['stop'] - km_data['start']

# Initialize Kaplan-Meier estimators
kmf_uv = KaplanMeierFitter()
kmf_vx = KaplanMeierFitter()
fig = go.Figure()

# Plot KM curves for each group
for group, label, color in zip(['Unvaccinated', 'Vaccinated'], ['Unvaccinated', 'Vaccinated'], ['blue', 'red']):
    mask = km_data['group'] == group
    kmf = KaplanMeierFitter()
    kmf.fit(durations=km_data.loc[mask, 'duration'], event_observed=km_data.loc[mask, 'event'], label=label)
    
    fig.add_trace(go.Scatter(
        x=kmf.survival_function_.index,
        y=kmf.survival_function_[label],
        mode='lines',
        name=label,
        line=dict(color=color)
    ))

# Update plot layout
fig.update_layout(
    title="Kaplan-Meier Survival Curve by Vaccination Status",
    xaxis_title="Time (days)",
    yaxis_title="Survival Probability",
    template='plotly_white',
    hovermode='x unified'
)

# Save KM plot as HTML
km_plot_path = OUTPUT_HTML.replace('.html', '_KM_survival.html')
fig.write_html(km_plot_path)
print(f"KM survival plot saved to {km_plot_path}")

# === Cleanup ===

# Restore original stdout/stderr and close log file
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()
