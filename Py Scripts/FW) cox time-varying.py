import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lifelines import CoxTimeVaryingFitter, KaplanMeierFitter
from scipy.integrate import simps  # for numerical integration
import sys

"""
Time-Varying Cox Regression and Survival Analysis on Vaccination and Death Data

This script performs a Cox time-varying survival analysis on individual-level vaccination and death data.
It uses real or simulated vaccination timing to define dynamic exposure windows, and computes hazard ratios
and life-years saved between vaccinated and unvaccinated groups.

Major steps include:
1. Loading individual-level data including vaccination and death dates.
2. Preprocessing to calculate age and convert date fields to numeric day indices.
3. Creating time-varying segments for each individual with vaccinated/unvaccinated periods.
4. Fitting a Cox time-varying model using lifelines.
5. Plotting Kaplan-Meier survival curves for both exposure groups using Plotly.
6. Estimating life-years saved using numerical integration.

Required:
- Input CSV with 'Rok_narozeni', 'DatumUmrti', and 'Datum_1' to 'Datum_7'
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lifelines import CoxTimeVaryingFitter, KaplanMeierFitter
from scipy.integrate import simps  # for numerical integration
import sys

# === Constants for I/O and analysis configuration ===

# Select one of the input/output paths depending on the dataset being analyzed
# Uncomment this block for simulated data

#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 cox time-varying.html"
#OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70 cox time-varying.TXT"

# real czech data 
INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW) Vesely_106_202403141131_AG70 cox time-varying.html"
OUTPUT_TXT = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FW) cox time-varying\FW) Vesely_106_202403141131_AG70 cox time-varying.TXT"

START_DATE = pd.Timestamp('2020-01-01')  # Reference date for calculating relative day indices
REFERENCE_YEAR = 2023                   # Used to calculate age from year of birth
MAX_AGE = 113                           # Age filtering threshold
LAG_DAYS = 0                            # Immunization lag (e.g., 14 days) after vaccination

original_stdout = sys.stdout  # Backup original stdout

# === Logging Setup ===

class Tee:
    """Helper class to duplicate stdout/stderr to both console and log file"""
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

# === Load CSV and preprocess ===

# Define vaccination date columns
dose_cols = [f'Datum_{i}' for i in range(1, 8)]

# Load only necessary columns and parse dates
df = pd.read_csv(
    INPUT_CSV,
    usecols=['Rok_narozeni', 'DatumUmrti'] + dose_cols,
    parse_dates=['DatumUmrti'] + dose_cols,
    dayfirst=False
)

# Normalize column names
df.columns = [col.lower().strip() for col in df.columns]
dose_cols = [col.lower() for col in dose_cols]

# Calculate age and drop invalid rows
df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']
df = df[df['age'].between(0, MAX_AGE)].copy()
df.dropna(subset=['age'], inplace=True)

# Override: If the dataset is already for age 70 group, force all rows to age 70
df['age'] = 70

# === Convert all date columns to integer day numbers ===
to_day = lambda col: (df[col] - START_DATE).dt.days
df['death_day'] = to_day('datumumrti')
for col in dose_cols:
    df[col + '_day'] = to_day(col)

# Identify first dose day and whether any dose exists
dose_day_cols = [col + '_day' for col in dose_cols]
df['first_dose_day'] = df[dose_day_cols].min(axis=1, skipna=True)
df['has_any_dose'] = df[dose_day_cols].notna().any(axis=1)

# Define the last measurement day (for censoring)
END_MEASURE = int(df['death_day'].dropna().max())
df['end_day'] = df['death_day'].fillna(END_MEASURE)
print(f"END_MEASURE (max death day): {END_MEASURE}")

# === Create Time-Varying Format for Cox Model ===

tv_data = []
for idx, row in df.iterrows():
    pid = idx  # Use index as patient ID
    death_day = row['death_day']
    end_day = row['end_day']
    dose_day = row['first_dose_day']
    has_dose = row['has_any_dose']

    # Segment 1: Pre-vaccination period (until dose+lag or death)
    unvax_stop = min(end_day, dose_day + LAG_DAYS) if has_dose else end_day
    tv_data.append({
        'id': pid,
        'start': 0,
        'stop': unvax_stop,
        'event': int(death_day == unvax_stop),
        'vaccinated': 0,
        't': unvax_stop
    })

    # Segment 2: Post-vaccination period (dose+lag to death or censoring)
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

# Adjust rows where start == stop and event==1 to avoid 0-duration intervals
tv_df.loc[(tv_df["start"] == tv_df["stop"]) & (tv_df["event"] == 1), "stop"] += 0.5

# Add a time-dependent interaction term
tv_df['vaccinated_time'] = tv_df['vaccinated'] * (tv_df['stop'] - tv_df['start'])

# Clean NaNs and infinities before fitting
tv_df.replace([np.inf, -np.inf], np.nan, inplace=True)
tv_df.dropna(inplace=True)

# === Fit Cox Time-Varying Model ===

ctv = CoxTimeVaryingFitter(penalizer=0.1)
ctv.fit(tv_df, id_col="id", start_col="start", stop_col="stop", event_col="event", show_progress=True)

# Print model summary
print(ctv.summary)

# === Compute Hazard Ratios with Confidence Intervals ===

hr = ctv.hazard_ratios_
ci = ctv.confidence_intervals_

# Identify columns for confidence interval bounds
lower_col = [col for col in ci.columns if "lower" in col.lower()][0]
upper_col = [col for col in ci.columns if "upper" in col.lower()][0]

# Output HRs with 95% CI
print("\nHazard Ratios and 95% Confidence Intervals:")
for cov in hr.index:
    lb = ci.loc[cov, lower_col]
    ub = ci.loc[cov, upper_col]
    print(f"{cov}: HR = {hr[cov]:.3f} (95% CI: {np.exp(lb):.3f} - {np.exp(ub):.3f})")

# === Plot Kaplan-Meier Survival Curves using Plotly ===

# Fit KM model to unvaccinated intervals
kmf_uvx = KaplanMeierFitter()
mask_uvx = tv_df['vaccinated'] == 0
durations_uvx = tv_df.loc[mask_uvx, 'stop'] - tv_df.loc[mask_uvx, 'start']
events_uvx = tv_df.loc[mask_uvx, 'event']
kmf_uvx.fit(durations=durations_uvx, event_observed=events_uvx, label="Unvaccinated")

# Fit KM model to vaccinated intervals
kmf_vx = KaplanMeierFitter()
mask_vx = tv_df['vaccinated'] == 1
durations_vx = tv_df.loc[mask_vx, 'stop'] - tv_df.loc[mask_vx, 'start']
events_vx = tv_df.loc[mask_vx, 'event']
kmf_vx.fit(durations=durations_vx, event_observed=events_vx, label="Vaccinated")

fig = go.Figure()

# Add unvaccinated survival trace
fig.add_trace(go.Scatter(
    x=kmf_uvx.survival_function_.index,
    y=kmf_uvx.survival_function_['Unvaccinated'],
    mode='lines',
    name='Unvaccinated',
    line=dict(color='blue')
))

# Add vaccinated survival trace
fig.add_trace(go.Scatter(
    x=kmf_vx.survival_function_.index,
    y=kmf_vx.survival_function_['Vaccinated'],
    mode='lines',
    name='Vaccinated',
    line=dict(color='red')
))

# Layout adjustments
fig.update_layout(
    title="Survival curves by vaccination state (Time-Varying Cox Model)",
    xaxis_title="Days under exposure",
    yaxis_title="Survival probability",
    template='plotly_white',
    hovermode="x unified"
)

# === Calculate Life Years Saved by Integration ===

# Define max integration limit across both curves
max_day_uvx = kmf_uvx.survival_function_.index.max()
max_day_vx = kmf_vx.survival_function_.index.max()
max_day = min(END_MEASURE, max_day_uvx, max_day_vx)

# Truncate survival curves to same time range
surv_uvx = kmf_uvx.survival_function_.loc[:max_day, 'Unvaccinated']
surv_vx = kmf_vx.survival_function_.loc[:max_day, 'Vaccinated']

# Time axis
time_uvx = surv_uvx.index.values
time_vx = surv_vx.index.values

# Numerical integration to compute expected survival time
expected_surv_uvx = simps(surv_uvx.values, time_uvx)
expected_surv_vx = simps(surv_vx.values, time_vx)

# Difference in expected survival time = life-years saved
life_years_saved = (expected_surv_vx - expected_surv_uvx) / 365
print(f"Life years saved (vaccinated vs unvaccinated) up to day {max_day}: {life_years_saved:.4f} years")

# Add annotation to survival plot
fig.add_annotation(
    x=max_day * 0.7,
    y=0.1,
    text=f"Life Years Saved: {life_years_saved:.3f} years",
    showarrow=False,
    font=dict(size=14, color="green"),
    bgcolor="rgba(255,255,255,0.8)"
)

# Save final plot as HTML
fig.write_html(OUTPUT_HTML)
print(f"Plot saved to {OUTPUT_HTML} with life years saved annotation")

# === Cleanup: Restore original stdout/stderr and close log ===
sys.stdout = original_stdout
sys.stderr = original_stderr
log_file.close()
