import pandas as pd
import numpy as np
from lifelines import CoxTimeVaryingFitter, KaplanMeierFitter
import plotly.graph_objects as go
import sys

# === Constants ===

#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FX) TTE per dose\FX-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 TTE.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FX) TTE per dose\FX-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 TTE.TXT"

INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FX) TTE per dose\FX) Vesely_106_202403141131_AG70 TTE.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FX) TTE per dose\FX) Vesely_106_202403141131_AG70 TTE.TXT"


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

# === Define end of observation ===
END_MEASURE = int(df['death_day'].dropna().max())
print(f"END_MEASURE (max death day): {END_MEASURE}")
df['end_day'] = df['death_day'].fillna(END_MEASURE)
df['event'] = (~df['death_day'].isna()).astype(int)

# === Prepare time-varying exposure records ===
records = []

for _, row in df.iterrows():
    pid = row.name
    death_day = row['death_day'] if not np.isnan(row['death_day']) else np.inf
    end_day = int(row['end_day'])

    dose_days = []
    for i, col in enumerate(DOSE_COLS):
        day = row.get(col + '_day')
        if not np.isnan(day):
            dose_days.append((int(day) + IMMUNITY_LAG, i + 1))

    dose_days.sort()
    current_day = 0
    current_dose = 0

    for start_day, dose_number in dose_days:
        if start_day >= end_day:
            break

        stop_day = min(start_day, end_day)
        event = int(death_day == stop_day)

        if stop_day > current_day:
            records.append({
                'id': pid,
                'start': current_day,
                'stop': stop_day,
                'dose_number': current_dose,
                'event': event
            })

        current_day = start_day
        current_dose = dose_number

        if death_day <= current_day:
            break

    if current_day < end_day:
        event = int(death_day == end_day)
        records.append({
            'id': pid,
            'start': current_day,
            'stop': end_day,
            'dose_number': current_dose,
            'event': event
        })

tte_df = pd.DataFrame(records)

# Fix zero-length intervals with event
tte_df.loc[
    (tte_df["start"] == tte_df["stop"]) & (tte_df["event"] == 1),
    "stop"
] += 0.5

# === Dummy coding for doses (baseline = dose 0) ===
tte_df = pd.get_dummies(tte_df, columns=["dose_number"], prefix="dose", drop_first=True)

# === Fit Cox Time-Varying Model ===
ctv = CoxTimeVaryingFitter(penalizer=0.1)
ctv.fit(tte_df, id_col="id", start_col="start", stop_col="stop", event_col="event")
ctv.print_summary()

# === Plot Survival Curves by Final Dose with Plotly ===
kmf = KaplanMeierFitter()
fig = go.Figure()

for dose in range(0, 8):
    if dose == 0:
        group = tte_df[
            ~tte_df[[col for col in tte_df.columns if col.startswith("dose_")]].any(axis=1)
        ]
    else:
        dose_col = f'dose_{dose}'
        if dose_col not in tte_df.columns:
            continue
        group = tte_df[
            (tte_df[dose_col] == 1) &
            (tte_df[[col for col in tte_df.columns if col.startswith("dose_") and col != dose_col]].sum(axis=1) == 0)
        ]
    
    if group.empty:
        continue

    durations = group["stop"] - group["start"]
    events = group["event"]
    label = f"Dose {dose}"

    kmf.fit(durations=durations, event_observed=events, label=label)
    survival_df = kmf.survival_function_.reset_index()
    
    fig.add_trace(go.Scatter(
        x=survival_df["timeline"],
        y=survival_df[kmf._label],
        mode='lines',
        name=label
    ))

fig.update_layout(
    title="Survival Curves Stratified by Final Dose (Kaplan-Meier)",
    xaxis_title="Days since start",
    yaxis_title="Survival Probability",
    legend_title="Dose Number",
    template="plotly_white"
)

# === Save to HTML ===
fig.write_html(OUTPUT_HTML)
print(f"Plot saved to: {OUTPUT_HTML}")

# close logging console and restore original streams at end
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()