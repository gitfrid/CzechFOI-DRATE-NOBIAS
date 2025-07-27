import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lifelines import CoxTimeVaryingFitter, KaplanMeierFitter
from scipy.integrate import simps  # for numerical integration
import sys

# === Constants ===
#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 cox time-varying.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 cox time-varying.TXT"

INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW) Vesely_106_202403141131_AG70 cox time-varying.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW) Vesely_106_202403141131_AG70 cox time-varying.TXT"

START_DATE = pd.Timestamp('2020-01-01')
REFERENCE_YEAR = 2023
MAX_AGE = 113
LAG_DAYS = 0  # Immunization starts 14 days after vaccination
original_stdout = sys.stdout # Save original stdout

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
    dose_day = row['first_dose_day']
    has_dose = row['has_any_dose']

    # Segment 1: Unvaccinated period (until immunization or death)
    unvax_stop = min(end_day, dose_day + LAG_DAYS) if has_dose else end_day
    tv_data.append({
        'id': pid,
        'start': 0,
        'stop': unvax_stop,
        'event': int(death_day == unvax_stop),
        'vaccinated': 0,
        't': unvax_stop
    })

    # Segment 2: Vaccinated period (from immunization to death or censoring)
    if has_dose and dose_day + LAG_DAYS < end_day:
        tv_data.append({
            'id': pid,
            'start': dose_day + LAG_DAYS,
            'stop': end_day,
            'event': int(death_day == end_day),
            'vaccinated': 1,
            't': end_day
        })

tv_df = pd.DataFrame(tv_data)

# Correct intervals where start == stop and event==1 by extending stop a bit
tv_df.loc[(tv_df["start"] == tv_df["stop"]) & (tv_df["event"] == 1), "stop"] += 0.5

# Add time-dependent interaction term vaccinated_time = vaccinated * (stop - start)
tv_df['vaccinated_time'] = tv_df['vaccinated'] * (tv_df['stop'] - tv_df['start'])

# Remove NaNs and infinities
tv_df.replace([np.inf, -np.inf], np.nan, inplace=True)
tv_df.dropna(inplace=True)

# === Fit Cox Time-Varying Model ===
ctv = CoxTimeVaryingFitter(penalizer=0.1)
ctv.fit(tv_df, id_col="id", start_col="start", stop_col="stop", event_col="event", show_progress=True)

print(ctv.summary)

# === Hazard Ratios and Confidence Intervals (Corrected) ===
hr = ctv.hazard_ratios_
ci = ctv.confidence_intervals_  # Has "lower-bound" and "upper-bound"

lower_col = [col for col in ci.columns if "lower" in col.lower()][0]
upper_col = [col for col in ci.columns if "upper" in col.lower()][0]

print("\nHazard Ratios and 95% Confidence Intervals:")
for cov in hr.index:
    lb = ci.loc[cov, lower_col]
    ub = ci.loc[cov, upper_col]
    print(f"{cov}: HR = {hr[cov]:.3f} (95% CI: {np.exp(lb):.3f} - {np.exp(ub):.3f})")

# === Plot Kaplan-Meier survival curves by vaccination status using Plotly ===
kmf_vx = KaplanMeierFitter()
kmf_uvx = KaplanMeierFitter()

mask_uvx = tv_df['vaccinated'] == 0
durations_uvx = tv_df.loc[mask_uvx, 'stop'] - tv_df.loc[mask_uvx, 'start']
events_uvx = tv_df.loc[mask_uvx, 'event']
kmf_uvx.fit(durations=durations_uvx, event_observed=events_uvx, label="Unvaccinated")

mask_vx = tv_df['vaccinated'] == 1
durations_vx = tv_df.loc[mask_vx, 'stop'] - tv_df.loc[mask_vx, 'start']
events_vx = tv_df.loc[mask_vx, 'event']
kmf_vx.fit(durations=durations_vx, event_observed=events_vx, label="Vaccinated")

fig = go.Figure()

# Add unvaccinated survival function
fig.add_trace(go.Scatter(
    x=kmf_uvx.survival_function_.index,
    y=kmf_uvx.survival_function_['Unvaccinated'],
    mode='lines',
    name='Unvaccinated',
    line=dict(color='blue')
))

# Add vaccinated survival function
fig.add_trace(go.Scatter(
    x=kmf_vx.survival_function_.index,
    y=kmf_vx.survival_function_['Vaccinated'],
    mode='lines',
    name='Vaccinated',
    line=dict(color='red')
))

fig.update_layout(
    title="Survival curves by vaccination state (Time-Varying Cox Model)",
    xaxis_title="Days under exposure",
    yaxis_title="Survival probability",
    template='plotly_white',
    hovermode="x unified"
)

# === Calculate Life Years Saved ===
max_day_uvx = kmf_uvx.survival_function_.index.max()
max_day_vx = kmf_vx.survival_function_.index.max()
max_day = min(END_MEASURE, max_day_uvx, max_day_vx)

surv_uvx = kmf_uvx.survival_function_.loc[:max_day, 'Unvaccinated']
surv_vx = kmf_vx.survival_function_.loc[:max_day, 'Vaccinated']

time_uvx = surv_uvx.index.values
time_vx = surv_vx.index.values

expected_surv_uvx = simps(surv_uvx.values, time_uvx)
expected_surv_vx = simps(surv_vx.values, time_vx)

life_years_saved = (expected_surv_vx - expected_surv_uvx) / 365

print(f"Life years saved (vaccinated vs unvaccinated) up to day {max_day}: {life_years_saved:.4f} years")

# Add annotation with Life Years Saved on the plot
fig.add_annotation(
    x=max_day * 0.7,
    y=0.1,
    text=f"Life Years Saved: {life_years_saved:.3f} years",
    showarrow=False,
    font=dict(size=14, color="green"),
    bgcolor="rgba(255,255,255,0.8)"
)

fig.write_html(OUTPUT_HTML)
print(f"Plot saved to {OUTPUT_HTML} with life years saved annotation")

# close logging console and restore original streams at end
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()
