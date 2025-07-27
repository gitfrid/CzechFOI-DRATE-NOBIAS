import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lifelines import CoxTimeVaryingFitter, KaplanMeierFitter
from scipy.integrate import simps  # for numerical integration
import sys

# === Constants ===

# INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
# OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FY) cox time-varying per Dose\FY-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 cox time-varying per dose.html"
# OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FY) cox time-varying per Dose\FY-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 cox time-varying per dose.TXT"


INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FY) cox time-varying per Dose\FY) real data Vesely_106_202403141131_AG70 cox time-varying per dose.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FY) cox time-varying per Dose\FY) real data Vesely_106_202403141131_AG70 cox time-varying per dose.TXT"

START_DATE = pd.Timestamp('2020-01-01')
REFERENCE_YEAR = 2023
MAX_AGE = 113
LAG_DAYS = 0  # Immunization starts 14 days after vaccination

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

# === Load CSV ===
dose_cols = [f'Datum_{i}' for i in range(1, 8)]
df = pd.read_csv(
    INPUT_CSV,
    usecols=['Rok_narozeni', 'DatumUmrti'] + dose_cols,
    parse_dates=['DatumUmrti'] + dose_cols,
    dayfirst=False
)
df.columns = [col.lower().strip() for col in df.columns]
dose_cols = [col.lower() for col in dose_cols]

# === Calculate age and filter ===
df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']
df = df[df['age'].between(0, MAX_AGE)].copy()
df.dropna(subset=['age'], inplace=True)

# Optional: fix age to 70 if dataset is only for age 70 group
df['age'] = 70

# === Convert dates to day numbers ===
to_day = lambda col: (df[col] - START_DATE).dt.days
df['death_day'] = to_day('datumumrti')
for col in dose_cols:
    df[col + '_day'] = to_day(col)

dose_day_cols = [col + '_day' for col in dose_cols]
df['first_dose_day'] = df[dose_day_cols].min(axis=1, skipna=True)
df['has_any_dose'] = df[dose_day_cols].notna().any(axis=1)

END_MEASURE = int(df['death_day'].dropna().max())
df['end_day'] = df['death_day'].fillna(END_MEASURE)

print(f"END_MEASURE (max death day): {END_MEASURE}")

# === Prepare time-varying data for Cox model ===
tv_data = []
for idx, row in df.iterrows():
    pid = idx
    death_day = row['death_day']
    end_day = row['end_day']

    # For stratification by dose number, we create intervals per dose:
    dose_days = [row[col + '_day'] for col in dose_cols]
    # Replace NaN with inf for sorting
    dose_days_sorted = sorted([d if not pd.isna(d) else np.inf for d in dose_days])

    last_start = 0
    for dose_num, dose_day in enumerate(dose_days_sorted):
        # If dose_day is inf, means no more doses
        if dose_day == np.inf or dose_day + LAG_DAYS > end_day:
            break
        # Unvaccinated or previous dose interval
        tv_data.append({
            'id': pid,
            'start': last_start,
            'stop': dose_day + LAG_DAYS,
            'event': int(death_day == dose_day + LAG_DAYS),
            'dose_num': dose_num,  # dose_num = 0 for unvaccinated, 1 for first dose, etc.
            't': dose_day + LAG_DAYS
        })
        last_start = dose_day + LAG_DAYS

    # Last interval: from last dose or start to end_day
    if last_start < end_day:
        tv_data.append({
            'id': pid,
            'start': last_start,
            'stop': end_day,
            'event': int(death_day == end_day),
            'dose_num': len([d for d in dose_days_sorted if d != np.inf]),
            't': end_day
        })

tv_df = pd.DataFrame(tv_data)

# Correct intervals where start == stop and event==1 by extending stop a bit
tv_df.loc[(tv_df["start"] == tv_df["stop"]) & (tv_df["event"] == 1), "stop"] += 0.5

# Remove NaNs and infinities
tv_df.replace([np.inf, -np.inf], np.nan, inplace=True)
tv_df.dropna(inplace=True)

# === Fit Cox Time-Varying Model ===
ctv = CoxTimeVaryingFitter(penalizer=0.1)
ctv.fit(tv_df, id_col="id", start_col="start", stop_col="stop", event_col="event", show_progress=True)

print(ctv.summary)

# === Plot Kaplan-Meier survival curves stratified by dose number using Plotly ===
fig = go.Figure()
colors = ['black', 'blue', 'red', 'green', 'orange', 'purple', 'brown', 'cyan']  # Up to 8 dose states (including unvaccinated=0)
labels = ['Unvaccinated'] + [f'Dose {i}' for i in range(1, 8)]

for dose_num in sorted(tv_df['dose_num'].unique()):
    kmf = KaplanMeierFitter()
    mask = tv_df['dose_num'] == dose_num
    durations = tv_df.loc[mask, 'stop'] - tv_df.loc[mask, 'start']
    events = tv_df.loc[mask, 'event']
    label = labels[dose_num] if dose_num < len(labels) else f'Dose {dose_num}'
    kmf.fit(durations=durations, event_observed=events, label=label)

    fig.add_trace(go.Scatter(
        x=kmf.survival_function_.index,
        y=kmf.survival_function_[label],
        mode='lines',
        name=label,
        line=dict(color=colors[dose_num] if dose_num < len(colors) else None)
    ))

fig.update_layout(
    title="Survival Curves Stratified by Dose Number",
    xaxis_title="Days under exposure",
    yaxis_title="Survival probability",
    template='plotly_white',
    hovermode="x unified"
)

fig.write_html(OUTPUT_HTML)
print(f"Plot saved to {OUTPUT_HTML} with survival curves stratified by dose number")

# close logging console and restore original streams at end
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()
