import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import KaplanMeierFitter
import plotly.graph_objects as go
import sys

# === Constants and input ===
INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.TXT"


#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP) real data Vesely_106_202403141131_AG70 poisson.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FP) poisson speedup\FP) real data Vesely_106_202403141131_AG70 poisson.TXT"


START_DATE = pd.Timestamp('2020-01-01')
MAX_AGE = 113
REFERENCE_YEAR = 2023

# --- Logging helper to tee stdout and stderr ---
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

# Setup logging to console and file simultaneously
original_stdout = sys.stdout
original_stderr = sys.stderr
log_file = open(OUTPUT_TXT, "w", encoding="utf-8")
tee = Tee(sys.stdout, log_file)
sys.stdout = tee
sys.stderr = tee

# === Load and Prepare Data ===
dose_date_cols = [f'Datum_{i}' for i in range(1, 8)]
needed_cols = ['Rok_narozeni', 'DatumUmrti'] + dose_date_cols

print("Loading data...")
df = pd.read_csv(
    INPUT_CSV,
    usecols=needed_cols,
    parse_dates=['DatumUmrti'] + dose_date_cols,
    dayfirst=False,
    low_memory=False
)

# Lowercase column names
df.columns = [col.strip().lower() for col in df.columns]
dose_date_cols_lower = [col.lower() for col in dose_date_cols]

df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']
df = df[df['age'].between(0, MAX_AGE)].copy()

def to_day_number(date_series):
    return (date_series - START_DATE).dt.days

df['death_day'] = to_day_number(df['datumumrti'])
for col in dose_date_cols_lower:
    df[col + '_day'] = to_day_number(df[col])

# First dose day per individual
df['first_dose_day'] = df[[col + '_day' for col in dose_date_cols_lower]].min(axis=1, skipna=True)
df['has_any_dose'] = df[[col + '_day' for col in dose_date_cols_lower]].notna().any(axis=1)

END_MEASURE = int(df['death_day'].dropna().max())
print(f"END_MEASURE (max death day): {END_MEASURE}")

# Define end of follow-up day (death or censoring)
df['end_day'] = df['death_day'].fillna(END_MEASURE)

# === Expand each person into daily rows ===

def expand_person_days(row):
    """Create daily records for one person from day 0 to end_day inclusive.
    vaccinated = 1 if day >= first_dose_day and person has any dose, else 0."""
    end = int(row['end_day'])
    days = np.arange(end + 1)  # from 0 to end_day inclusive
    vaccinated = np.zeros_like(days, dtype=int)
    if pd.notna(row['first_dose_day']) and row['has_any_dose']:
        vaccinated[days >= row['first_dose_day']] = 1
    return pd.DataFrame({
        'age': row['age'],
        'day': days,
        'vaccinated': vaccinated,
        'death': 0
    })

print("Expanding person-day data... this may take a few minutes for large data.")

person_days_list = []
for i, row in df.iterrows():
    person_days_list.append(expand_person_days(row))

person_days = pd.concat(person_days_list, ignore_index=True)

# Mark deaths in person_days
death_df = df.loc[df['death_day'].notna(), ['age', 'death_day']]
for _, drow in death_df.iterrows():
    mask = (person_days['age'] == drow['age']) & (person_days['day'] == drow['death_day'])
    person_days.loc[mask, 'death'] = 1

# === Aggregate counts by age, day, vaccination status ===
print("Aggregating data by age, day, vaccination status...")
agg = person_days.groupby(['age', 'day', 'vaccinated']).agg(
    deaths=('death', 'sum'),
    person_days=('day', 'count')
).reset_index()

# Offset for Poisson regression = log(person_days)
agg['offset'] = np.log(agg['person_days'])

# Center age variable for modeling
agg['age_c'] = agg['age'] - agg['age'].mean()

# Design matrix with intercept, vaccination, centered age
X = sm.add_constant(agg[['vaccinated', 'age_c']])

# Fit Poisson regression model
print("Fitting Poisson regression model...")
model = sm.GLM(agg['deaths'], X, offset=agg['offset'], family=sm.families.Poisson())
result = model.fit()

print(result.summary())

params = result.params
conf = result.conf_int()
irr = np.exp(params)
irr_conf_lower = np.exp(conf[0])
irr_conf_upper = np.exp(conf[1])

print("\nIncidence Rate Ratios (IRRs):")
print(f"Intercept: {irr['const']:.3f}")
print(f"Vaccinated vs Unvaccinated: {irr['vaccinated']:.3f} "
      f"(95% CI: {irr_conf_lower['vaccinated']:.3f} - {irr_conf_upper['vaccinated']:.3f})")

# === Kaplan-Meier Survival Analysis ===

print("Preparing Kaplan-Meier survival data...")

# Unvaccinated period: from day 0 to first dose day or end_day
unvaccinated = df[['age', 'end_day', 'first_dose_day', 'death_day']].copy()
unvaccinated['start'] = 0
unvaccinated['stop'] = unvaccinated['first_dose_day'].fillna(unvaccinated['end_day'])
unvaccinated['event'] = ((unvaccinated['death_day'] <= unvaccinated['stop']) & (unvaccinated['death_day'].notna())).astype(int)
unvaccinated['group'] = 'Unvaccinated'
unvaccinated = unvaccinated[unvaccinated['stop'] > unvaccinated['start']]

# Vaccinated period: from first dose day to end_day
vaccinated = df[df['first_dose_day'].notna()][['age', 'end_day', 'first_dose_day', 'death_day']].copy()
vaccinated['start'] = vaccinated['first_dose_day']
vaccinated['stop'] = vaccinated['end_day']
vaccinated['event'] = ((vaccinated['death_day'] >= vaccinated['start']) & (vaccinated['death_day'].notna())).astype(int)
vaccinated['group'] = 'Vaccinated'
vaccinated = vaccinated[vaccinated['stop'] > vaccinated['start']]

km_data = pd.concat([unvaccinated, vaccinated], ignore_index=True)
km_data['duration'] = km_data['stop'] - km_data['start']

# Plot Kaplan-Meier survival curves by group using lifelines + plotly
fig = go.Figure()
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

fig.update_layout(
    title="Kaplan-Meier Survival Curve by Vaccination Status",
    xaxis_title="Time (days)",
    yaxis_title="Survival Probability",
    template='plotly_white',
    hovermode='x unified'
)

km_plot_path = OUTPUT_HTML.replace('.html', '_KM_survival.html')
fig.write_html(km_plot_path)
print(f"KM survival plot saved to {km_plot_path}")

# Close log file and restore stdout/stderr
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()
print("Script completed.")
