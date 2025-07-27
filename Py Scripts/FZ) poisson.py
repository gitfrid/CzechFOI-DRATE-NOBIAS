import pandas as pd
import numpy as np
import statsmodels.api as sm
from lifelines import KaplanMeierFitter
import plotly.graph_objects as go
import sys

# === Constants and input ===
INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 poisson.TXT"

#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ) real data Vesely_106_202403141131_AG70 poisson.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FZ) poisson\FZ) real data Vesely_106_202403141131_AG70 poisson.TXT"

START_DATE = pd.Timestamp('2020-01-01')
MAX_AGE = 113
REFERENCE_YEAR = 2023

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

# Save original stdout/stderr to restore later if needed
original_stdout = sys.stdout
original_stderr = sys.stderr

log_file_path = fr"{OUTPUT_TXT}"
log_file = open(log_file_path, "w", encoding="utf-8")

# Set up teeing
tee = Tee(sys.stdout, log_file)
sys.stdout = tee
sys.stderr = tee

# === Load and Prepare Data ===
dose_date_cols = [f'Datum_{i}' for i in range(1, 8)]
needed_cols = ['Rok_narozeni', 'DatumUmrti'] + dose_date_cols

df = pd.read_csv(
    INPUT_CSV,
    usecols=needed_cols,
    parse_dates=['DatumUmrti'] + dose_date_cols,
    dayfirst=False,
    low_memory=False
)

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

df['first_dose_day'] = df[[col + '_day' for col in dose_date_cols_lower]].min(axis=1, skipna=True)
df['has_any_dose'] = df[[col + '_day' for col in dose_date_cols_lower]].notna().any(axis=1)

# === Define observation period ===
END_MEASURE = int(df['death_day'].dropna().max())
print(f"END_MEASURE (max death day): {END_MEASURE}")

# For those who died, end_day = death_day; for alive, end_day = END_MEASURE (censored)
df['end_day'] = df['death_day'].fillna(END_MEASURE)

# === Create follow-up time and vaccination status with zero lag post-dose ===

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

person_days_list = []
for idx, row in df.iterrows():
    person_days_list.append(expand_person_days(row))

person_days = pd.concat(person_days_list, ignore_index=True)

death_rows = df.loc[df['death_day'].notna(), ['age', 'death_day']]
for _, death_row in death_rows.iterrows():
    mask = (person_days['age'] == death_row['age']) & (person_days['day'] == death_row['death_day'])
    person_days.loc[mask, 'death'] = 1

# Aggregate by age, day, vaccinated status
agg = person_days.groupby(['age', 'day', 'vaccinated']).agg(
    deaths=('death', 'sum'),
    person_days=('day', 'count')
).reset_index()

agg['offset'] = np.log(agg['person_days'])
agg['age_c'] = agg['age'] - agg['age'].mean()
X = sm.add_constant(agg[['vaccinated', 'age_c']])

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
print(f"Vaccinated vs Unvaccinated: {irr['vaccinated']:.3f} (95% CI: {irr_conf_lower['vaccinated']:.3f} - {irr_conf_upper['vaccinated']:.3f})")

# === Kaplan-Meier Survival Plot ===

# Unvaccinated period (from 0 to first dose or censor/death)
unvaccinated = df[['age', 'end_day', 'first_dose_day', 'death_day']].copy()
unvaccinated['start'] = 0
unvaccinated['stop'] = unvaccinated['first_dose_day'].fillna(unvaccinated['end_day'])
unvaccinated['event'] = (unvaccinated['death_day'] <= unvaccinated['stop']) & (unvaccinated['death_day'].notna())
unvaccinated['event'] = unvaccinated['event'].astype(int)
unvaccinated['group'] = 'Unvaccinated'
unvaccinated = unvaccinated[unvaccinated['stop'] > unvaccinated['start']]

# Vaccinated period (from first dose to end)
vaccinated = df[df['first_dose_day'].notna()][['age', 'end_day', 'first_dose_day', 'death_day']].copy()
vaccinated['start'] = vaccinated['first_dose_day']
vaccinated['stop'] = vaccinated['end_day']
vaccinated['event'] = (vaccinated['death_day'] >= vaccinated['start']) & (vaccinated['death_day'].notna())
vaccinated['event'] = vaccinated['event'].astype(int)
vaccinated['group'] = 'Vaccinated'
vaccinated = vaccinated[vaccinated['stop'] > vaccinated['start']]

km_data = pd.concat([unvaccinated, vaccinated], ignore_index=True)
km_data['duration'] = km_data['stop'] - km_data['start']

kmf_uv = KaplanMeierFitter()
kmf_vx = KaplanMeierFitter()

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

# close logging console and restore original streams at end
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()