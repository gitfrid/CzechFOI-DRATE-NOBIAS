import pandas as pd
import numpy as np
from lifelines import CoxTimeVaryingFitter, KaplanMeierFitter
import plotly.graph_objects as go
import sys

# === Constants ===

INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FS) TTE\FS-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 TTE.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FS) TTE\FS-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 TTE.TXT"

#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FS) TTE\FS) Vesely_106_202403141131_AG70 TTE.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FS) TTE\FS) Vesely_106_202403141131_AG70 TTE.TXT"

START_DATE = pd.Timestamp('2020-01-01')
REFERENCE_YEAR = 2023
MAX_AGE = 113
IMMUNITY_LAG = 0  # days after dose until immunity starts
DOSE_COLS = [f'Datum_{i}' for i in range(1, 8)]

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

# === Load data ===
needed_cols = ['Rok_narozeni', 'DatumUmrti'] + DOSE_COLS
df = pd.read_csv(INPUT_CSV, usecols=needed_cols, parse_dates=['DatumUmrti'] + DOSE_COLS)
df.columns = [col.strip().lower() for col in df.columns]
DOSE_COLS = [col.lower() for col in DOSE_COLS]

df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']
df = df[df['age'].between(0, MAX_AGE)].copy()

# === Convert dates to days since START_DATE ===
def to_day_number(date_series):
    return (date_series - START_DATE).dt.days

df['death_day'] = to_day_number(df['datumumrti'])
for col in DOSE_COLS:
    df[col + '_day'] = to_day_number(df[col])

df['first_dose_day'] = df[[col + '_day' for col in DOSE_COLS]].min(axis=1, skipna=True)
df['has_any_dose'] = df[[col + '_day' for col in DOSE_COLS]].notna().any(axis=1)

# Define end of observation: max death day or administrative censoring day
END_MEASURE = int(df['death_day'].dropna().max())
print(f"END_MEASURE (max death day): {END_MEASURE}")

df['end_day'] = df['death_day'].fillna(END_MEASURE)
df['event'] = (~df['death_day'].isna()).astype(int)

# === Prepare Target Trial Emulation (TTE) structure ===
records = []

for _, row in df.iterrows():
    pid = row.name
    end = int(row['end_day'])
    event = row['event']
    death_day = row['death_day'] if not np.isnan(row['death_day']) else np.inf
    dose_day = row['first_dose_day']
    
    if np.isnan(dose_day):
        # Never vaccinated: entire follow-up unvaccinated
        records.append({
            'id': pid,
            'start': 0,
            'stop': end,
            'vaccinated': 0,
            'event': int(death_day == end)
        })
    else:
        immune_start = int(dose_day + IMMUNITY_LAG)
        
        # Interval 1: unvaccinated time [0, immune_start)
        unvax_stop = min(immune_start, end)
        event_unvax = int(death_day == unvax_stop)
        records.append({
            'id': pid,
            'start': 0,
            'stop': unvax_stop,
            'vaccinated': 0,
            'event': event_unvax
        })

        if immune_start < end:
            # Interval 2: vaccinated time [immune_start, end]
            event_vax = int(death_day == end)
            records.append({
                'id': pid,
                'start': immune_start,
                'stop': end,
                'vaccinated': 1,
                'event': event_vax
            })

tte_df = pd.DataFrame(records)

# Fix for zero-length intervals with event
tte_df.loc[
    (tte_df["start"] == tte_df["stop"]) & (tte_df["event"] == 1),
    "stop"
] += 0.5

# === Fit time-dependent Cox model ===
ctv = CoxTimeVaryingFitter(penalizer=0.1)
ctv.fit(tte_df, id_col="id", start_col="start", stop_col="stop", event_col="event")
ctv.print_summary()

# === Plot stratified survival curves by dose using Kaplan-Meier estimators ===

kmf_vax = KaplanMeierFitter()
kmf_unvax = KaplanMeierFitter()

# Aggregate data per individual for KM plot
# We use the maximum stop time per individual in vaccinated and unvaccinated states, with event if event happened during that state

# For unvaccinated:
unvax_df = tte_df[tte_df['vaccinated'] == 0].groupby('id').agg({
    'stop': 'max',
    'event': 'max'
}).reset_index()

# For vaccinated:
vax_df = tte_df[tte_df['vaccinated'] == 1].groupby('id').agg({
    'stop': 'max',
    'event': 'max'
}).reset_index()

# Fit KM curves
kmf_unvax.fit(durations=unvax_df['stop'], event_observed=unvax_df['event'], label='Unvaccinated')
kmf_vax.fit(durations=vax_df['stop'], event_observed=vax_df['event'], label='Vaccinated')

# Prepare Plotly figure
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=kmf_unvax.survival_function_.index,
    y=kmf_unvax.survival_function_['Unvaccinated'],
    mode='lines',
    name='Unvaccinated',
    line=dict(color='red')
))

fig.add_trace(go.Scatter(
    x=kmf_vax.survival_function_.index,
    y=kmf_vax.survival_function_['Vaccinated'],
    mode='lines',
    name='Vaccinated',
    line=dict(color='green')
))

fig.update_layout(
    title="Stratified Survival Curves by Dose (Vaccinated vs Unvaccinated)",
    xaxis_title="Days since Start",
    yaxis_title="Survival Probability",
    yaxis=dict(range=[0, 1]),
    template="plotly_white"
)

# Save interactive plot to HTML
fig.write_html(OUTPUT_HTML)
print(f"Interactive survival curves plot saved to {OUTPUT_HTML}")

# close logging console and restore original streams at end
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()