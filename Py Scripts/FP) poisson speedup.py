import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import KaplanMeierFitter
import plotly.graph_objects as go
import sys


# =============================================================================
# Script: Poisson Regression & Survival Analysis of Vaccination Data (Age ≤ 113)
#
# Description:
#   - Loads individual-level czech-FOI vaccination and mortality data
#   - Expands data to person-day format
#   - Flags daily vaccination status and death events
#   - Performs Poisson regression to estimate effect of vaccination on death risk
#   - Computes Kaplan-Meier survival curves for vaccinated vs unvaccinated
#
# Output:
#   - Console + text file log
#   - HTML survival curve plot
#
# Date: 29.07.2025
# =============================================================================


# === Constants and input ===
INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.TXT"

# Uncomment below to use real dataset instead of simulation
#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP) real data Vesely_106_202403141131_AG70 poisson.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP) real data Vesely_106_202403141131_AG70 poisson.TXT"

START_DATE = pd.Timestamp('2020-01-01')  # Origin date for day number calculation
MAX_AGE = 113                            # Maximum allowed age for inclusion
REFERENCE_YEAR = 2023                    # Reference year to calculate age

# === Tee logging ===
# Redirect stdout and stderr to both console and file
class Tee:
    def __init__(self, *files): self.files = files
    def write(self, data): [f.write(data) or f.flush() for f in self.files]
    def flush(self): [f.flush() for f in self.files]

log_file = open(OUTPUT_TXT, "w", encoding="utf-8")
sys.stdout = sys.stderr = Tee(sys.stdout, log_file)

# === Load and Prepare Data ===
dose_date_cols = [f'Datum_{i}' for i in range(1, 8)]
print("Loading data...")
df = pd.read_csv(INPUT_CSV, usecols=['Rok_narozeni', 'DatumUmrti'] + dose_date_cols, parse_dates=['DatumUmrti'] + dose_date_cols, dayfirst=False)
df.columns = [c.lower() for c in df.columns]
dose_date_cols = [c.lower() for c in dose_date_cols]

# Compute age from birth year
df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']
df = df[df['age'].between(0, MAX_AGE)]

# Convert all dates to integer day numbers since START_DATE
def to_day_number(dates): return (dates - START_DATE).dt.days
df['death_day'] = to_day_number(df['datumumrti'])
for col in dose_date_cols:
    df[col + '_day'] = to_day_number(df[col])

# Determine first dose day and dose status per individual
dose_day_cols = [col + '_day' for col in dose_date_cols]
df['first_dose_day'] = df[dose_day_cols].min(axis=1, skipna=True)
df['has_any_dose'] = df[dose_day_cols].notna().any(axis=1)

# Define follow-up end as death day or max death day if censored
df['end_day'] = df['death_day'].fillna(df['death_day'].max())
END_MEASURE = int(df['end_day'].max())
print(f"END_MEASURE (max death day): {END_MEASURE}")

# === Expand to person-day format ===
print("Expanding person-day data...")

rows = []
# Loop over each person and create daily records up to end_day
for row in df.itertuples(index=False):
    end_day = int(row.end_day)
    days = np.arange(end_day + 1)
    vaccinated = (pd.notna(row.first_dose_day) and row.has_any_dose)
    vax_mask = days >= row.first_dose_day if vaccinated else np.zeros_like(days, dtype=bool)
    rows.append(pd.DataFrame({
        'age': row.age,
        'day': days,
        'vaccinated': vax_mask.astype(int),
        'death': 0
    }))

person_days = pd.concat(rows, ignore_index=True)

# === Mark deaths efficiently without merge ===
# Create lookup set of (age, death_day) to quickly mark death events
death_idx = set(zip(df.loc[df['death_day'].notna(), 'age'], df.loc[df['death_day'].notna(), 'death_day']))
person_days['death'] = [
    1 if (age, day) in death_idx else 0
    for age, day in zip(person_days['age'], person_days['day'])
]

# === Aggregate and Model ===
print("Aggregating data...")

# Group by age, day, vaccination status and count deaths and person-days
agg = person_days.groupby(['age', 'day', 'vaccinated'], sort=False).agg(
    deaths=('death', 'sum'),
    person_days=('death', 'size')
).reset_index()

# Add offset and centered age for Poisson regression
agg['offset'] = np.log(agg['person_days'])
agg['age_c'] = agg['age'] - agg['age'].mean()
X = sm.add_constant(agg[['vaccinated', 'age_c']])

# Fit Poisson GLM model with log offset
print("Fitting Poisson regression model...")
model = sm.GLM(agg['deaths'], X, offset=agg['offset'], family=sm.families.Poisson())
result = model.fit()
print(result.summary())

# Compute Incidence Rate Ratios (IRR) with 95% confidence intervals
params = result.params
conf = result.conf_int()
irr = np.exp(params)
irr_conf = np.exp(conf)

print("\nIncidence Rate Ratios (IRRs):")
print(f"Intercept: {irr['const']:.3f}")
print(f"Vaccinated vs Unvaccinated: {irr['vaccinated']:.3f} "
      f"(95% CI: {irr_conf.loc['vaccinated', 0]:.3f} - {irr_conf.loc['vaccinated', 1]:.3f})")

# === Kaplan-Meier Survival Analysis ===
print("Preparing Kaplan-Meier survival data...")

# Unvaccinated period: from day 0 until first dose or end of follow-up
unvaccinated = df.copy()
unvaccinated['start'] = 0
unvaccinated['stop'] = df['first_dose_day'].fillna(df['end_day'])
unvaccinated['event'] = ((df['death_day'] <= unvaccinated['stop']) & df['death_day'].notna()).astype(int)
unvaccinated['group'] = 'Unvaccinated'
unvaccinated = unvaccinated[unvaccinated['stop'] > unvaccinated['start']]

# Vaccinated period: from first dose until end of follow-up
vaccinated = df[df['first_dose_day'].notna()].copy()
vaccinated['start'] = vaccinated['first_dose_day']
vaccinated['stop'] = vaccinated['end_day']
vaccinated['event'] = ((vaccinated['death_day'] >= vaccinated['start']) & vaccinated['death_day'].notna()).astype(int)
vaccinated['group'] = 'Vaccinated'
vaccinated = vaccinated[vaccinated['stop'] > vaccinated['start']]

# Combine both groups for survival analysis
km_data = pd.concat([unvaccinated, vaccinated], ignore_index=True)
km_data['duration'] = km_data['stop'] - km_data['start']

# === Plot Kaplan-Meier Curves ===
fig = go.Figure()
for group, color in zip(['Unvaccinated', 'Vaccinated'], ['blue', 'red']):
    mask = km_data['group'] == group
    kmf = KaplanMeierFitter()
    kmf.fit(km_data.loc[mask, 'duration'], km_data.loc[mask, 'event'], label=group)
    fig.add_trace(go.Scatter(
        x=kmf.survival_function_.index,
        y=kmf.survival_function_[group],
        mode='lines',
        name=group,
        line=dict(color=color)
    ))

fig.update_layout(
    title="Kaplan-Meier Survival Curve by Vaccination Status",
    xaxis_title="Time (days)",
    yaxis_title="Survival Probability",
    template='plotly_white',
    hovermode='x unified'
)

# Save survival curve as HTML
km_plot_path = OUTPUT_HTML.replace('.html', '_KM_survival.html')
fig.write_html(km_plot_path)
print(f"KM survival plot saved to {km_plot_path}")

# === Cleanup ===
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__
log_file.close()
print("Script completed.")
