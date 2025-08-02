import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lifelines import CoxTimeVaryingFitter
from scipy.integrate import simps
import sys

# === User input configuration ===
#INPUT_CSV = r"C:\github\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
#OUTPUT_HTML = r"C:\github\CzechFOI-DRATE-NOBIAS\Plot Results\FL) cox time-varying survival curves daily\FL-FG) case3_sim_deaths_sim_real_doses_with_constraint cox time-varying.html"
#OUTPUT_TXT = r"C:\github\CzechFOI-DRATE-NOBIAS\Plot Results\FL) cox time-varying survival curves daily\FL-FG) case3_sim_deaths_sim_real_doses_with_constraint cox time-varying.TXT"

INPUT_CSV = r"C:\github\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\github\CzechFOI-DRATE-NOBIAS\Plot Results\FL) cox time-varying survival curves daily\FL) Vesely_106_202403141131_AG70 cox time-varying.html"
OUTPUT_TXT = r"C:\github\CzechFOI-DRATE-NOBIAS\Plot Results\FL) cox time-varying survival curves daily\FL) Vesely_106_202403141131_AG70 cox time-varying.TXT"

"""
Cox Time-Varying Survival Analysis with Life Years Saved Estimation and Plotting
================================================================================

This script reads a vaccination and death dataset, preprocesses it to compute time-varying survival intervals
based on vaccine administration timing. It fits a Cox proportional hazards model using `lifelines`,
estimates hazard ratios, and computes separate survival curves for vaccinated and unvaccinated individuals.

It calculates the life years saved (area between the curves), and plots the predicted survival functions
using Plotly.

Assumptions:
------------
- Input CSV must contain birth year, death date, and up to 7 dose date columns it uses real world Czech FOI data
- Vaccination lag is configurable with `LAG_DAYS`.
- All individuals are forced to age 70 for analysis purposes.

Author: drifting - on a sea of forgotten teardrops
"""

# === USER CONFIGURATION ===
START_DATE = pd.Timestamp('2020-01-01')  # Origin date to count days from
REFERENCE_YEAR = 2023                    # Used for age calculation
MAX_AGE = 113                            # Max allowed age
LAG_DAYS = 0                             # Immunization lag in days
COX_PENALIZER = 3                        # Chosen based on best alignment with expected HR~1 in simulation; helps stabilize estimates and prevent overfitting/underfitting
AGE = 70                                 # Filter to certain AG for faster testing

# Columns for vaccine dose dates
dose_cols = [f'Datum_{i}' for i in range(1, 8)]

# === Logging both to console and text file ===
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

log_file = open(OUTPUT_TXT, "w", encoding="utf-8")
sys.stdout = Tee(sys.stdout, log_file)
sys.stderr = Tee(sys.stderr, log_file)

print("=== Starting Cox Time-Varying Survival Analysis ===")

# === Load dataset ===
print("Loading CSV...")
df = pd.read_csv(INPUT_CSV,
                 usecols=['Rok_narozeni', 'DatumUmrti'] + dose_cols,
                 parse_dates=['DatumUmrti'] + dose_cols,
                 dayfirst=False)

# Normalize column names
df.columns = [c.lower() for c in df.columns]
dose_cols = [col.lower() for col in dose_cols]

# Calculate age
df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
df['age'] = REFERENCE_YEAR - df['birth_year']
df = df[df['age'].between(0, MAX_AGE)].copy()
df.dropna(subset=['age'], inplace=True)

# Filter for AG 
df = df[df['age'] == AGE].copy()

# Convert all date columns to integer "days since START_DATE"
def to_day(col):
    return (df[col] - START_DATE).dt.days

df['death_day'] = to_day('datumumrti')
for col in dose_cols:
    df[col + '_day'] = to_day(col)

dose_day_cols = [col + '_day' for col in dose_cols]
df['first_dose_day'] = df[dose_day_cols].min(axis=1, skipna=True)
df['has_any_dose'] = df[dose_day_cols].notna().any(axis=1)

# Set the observation end day
END_MEASURE = int(df['death_day'].dropna().max())
df['end_day'] = df['death_day'].fillna(END_MEASURE)

print(f"Dataset filtered to age {MAX_AGE}, END_MEASURE (max death day): {END_MEASURE}")
print(f"Total rows after filtering: {df.shape[0]}")

# === Construct time-varying survival intervals ===
print("Preparing time-varying data for Cox model...")

tv_data = []
for idx, row in df.iterrows():
    pid = idx
    death_day = row['death_day']
    end_day = row['end_day']
    dose_day = row['first_dose_day']
    has_dose = row['has_any_dose']

    unvax_stop = min(end_day, dose_day + LAG_DAYS) if has_dose else end_day

    # Add unvaccinated period
    tv_data.append({
        'id': pid,
        'start': 0,
        'stop': unvax_stop,
        'event': int(death_day == unvax_stop),
        'vaccinated': 0
    })

    # Add vaccinated period if dose was given and subject survived lag
    if has_dose and dose_day + LAG_DAYS < end_day:
        tv_data.append({
            'id': pid,
            'start': dose_day + LAG_DAYS,
            'stop': end_day,
            'event': int(death_day == end_day),
            'vaccinated': 1
        })

tv_df = pd.DataFrame(tv_data)

# Fix zero-length intervals with event=1 by extending duration slightly
mask_zero_event = (tv_df['start'] == tv_df['stop']) & (tv_df['event'] == 1)
tv_df.loc[mask_zero_event, 'stop'] += 0.5

# Create a time-duration column for vaccinated segments (optional)
tv_df['vaccinated_time'] = tv_df['vaccinated'] * (tv_df['stop'] - tv_df['start'])

# Remove NaNs or infinities
tv_df.replace([np.inf, -np.inf], np.nan, inplace=True)
tv_df.dropna(inplace=True)

print(f"Prepared time-varying dataset with {tv_df.shape[0]} rows")

# === Fit the Cox model ===
print("Fitting Cox Time-Varying Model...")

ctv = CoxTimeVaryingFitter(penalizer=COX_PENALIZER)
ctv.fit(tv_df, id_col="id", start_col="start", stop_col="stop", event_col="event", show_progress=True)

print("\n=== Cox Model Summary ===")
print(ctv.summary)

# Extract hazard ratios and confidence intervals
hr = ctv.hazard_ratios_
ci = ctv.confidence_intervals_

lower_col = [c for c in ci.columns if "lower" in c.lower()][0]
upper_col = [c for c in ci.columns if "upper" in c.lower()][0]

print("\nHazard Ratios and 95% Confidence Intervals:")
for cov in hr.index:
    lb = ci.loc[cov, lower_col]
    ub = ci.loc[cov, upper_col]
    print(f"  {cov}: HR = {hr[cov]:.3f} (95% CI: {np.exp(lb):.3f} - {np.exp(ub):.3f})")

# === Compute survival curves for both groups ===
print("Calculating predicted daily survival curves manually...")

base_cumhaz = ctv.baseline_cumulative_hazard_

days = np.arange(0, END_MEASURE + 1)
base_cumhaz_interp = np.interp(days, base_cumhaz.index.values, base_cumhaz.values.flatten())
base_surv = np.exp(-base_cumhaz_interp)  # Baseline survival (unstratified)

# Extract model coefficients
coef = ctv.params_.to_dict()

# Survival curve generator
def survival_curve(vaccinated_flag):
    hr_vacc = np.exp(coef.get('vaccinated', 0) * vaccinated_flag)
    return base_surv ** hr_vacc

# Predict survival for vaccinated and unvaccinated
s_uvx = survival_curve(0)
s_vx = survival_curve(1)

# === Compute life years saved as area between survival curves ===
ly_uvx = simps(s_uvx, days)
ly_vx = simps(s_vx, days)
life_years_saved = ly_vx - ly_uvx

# === Debugging output of death events by vaccination ===
print(f"Life years saved by vaccination: {life_years_saved / 365.25:.4f} years")
events = tv_df['event'].astype(bool)
print("Variance vaccinated when event=1:", tv_df.loc[events, 'vaccinated'].var())
print("Variance vaccinated when event=0:", tv_df.loc[~events, 'vaccinated'].var())
print("Vaccinated counts when event=1:\n", tv_df.loc[events, 'vaccinated'].value_counts())
print("Vaccinated counts when event=0:\n", tv_df.loc[~events, 'vaccinated'].value_counts())
print("Deaths before vaccination + lag:", tv_df[(tv_df['vaccinated'] == 0) & (tv_df['event'] == 1)].shape[0])
print("Deaths after vaccination + lag:", tv_df[(tv_df['vaccinated'] == 1) & (tv_df['event'] == 1)].shape[0])
print("Count rows in time-varying dataset", tv_df.groupby(['event', 'vaccinated']).size())

event_mask = tv_df['event'] == 1
vx_events = tv_df[event_mask & (tv_df['vaccinated'] == 1)]
uvx_events = tv_df[event_mask & (tv_df['vaccinated'] == 0)]
print(f"\n=== DEBUG: Death Events by Vaccination Status ===")
print(f"Total death events: {event_mask.sum()}")
print(f"Vaccinated death events: {len(vx_events)}")
print(f"Unvaccinated death events: {len(uvx_events)}")
print("==============================================\n")

# === Plot survival curves using Plotly ===
print("Generating Plotly survival curve plot...")

fig = go.Figure()
fig.add_trace(go.Scatter(
    x=days,
    y=s_uvx,
    mode='lines',
    name='Unvaccinated',
    line=dict(color='red', width=3),
    hovertemplate='Day %{x}<br>Survival %{y:.4f}<extra>Unvaccinated</extra>'
))
fig.add_trace(go.Scatter(
    x=days,
    y=s_vx,
    mode='lines',
    name='Vaccinated',
    line=dict(color='green', width=3),
    hovertemplate='Day %{x}<br>Survival %{y:.4f}<extra>Vaccinated</extra>'
))

fig.update_layout(
    title=f"Cox Time-Varying Survival Curves (Age 70)<br>Life Years Saved: {life_years_saved:.2f} days",
    xaxis_title="Days since start",
    yaxis_title="Survival Probability",
    template="plotly_white"
)

fig.write_html(OUTPUT_HTML)
print(f"Plot saved to: {OUTPUT_HTML}")

print("=== Analysis Complete ===")

# === Cleanup logging ===
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__
log_file.close()
