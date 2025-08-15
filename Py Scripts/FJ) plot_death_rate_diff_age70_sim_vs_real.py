import pandas as pd
import numpy as np
import plotly.graph_objects as go
from lifelines import KaplanMeierFitter
from scipy.ndimage import gaussian_filter1d

# Input files
SIM_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
REAL_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FJ) bias vs observed vs adjusted KM death rate\AG70_bias_vs_observed_vs_adjusted_KM_death_rate.html"
OUTPUT_SURV_DIFF_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\FJ) bias vs observed vs adjusted KM death rate\AG70_KM_survival_difference.html"

START_DATE = pd.Timestamp('2020-01-01')
REFERENCE_YEAR = 2023
MAX_AGE = 113
LAG_DAYS = 0

dose_cols = [f'Datum_{i}' for i in range(1, 8)]

def preprocess_data(filepath):
    df = pd.read_csv(
        filepath,
        usecols=['Rok_narozeni', 'DatumUmrti'] + dose_cols,
        parse_dates=['DatumUmrti'] + dose_cols,
        dayfirst=False
    )
    df.columns = [col.lower().strip() for col in df.columns]
    dose_cols_lower = [col.lower() for col in dose_cols]

    df['birth_year'] = pd.to_numeric(df['rok_narozeni'], errors='coerce')
    df['age'] = REFERENCE_YEAR - df['birth_year']
    df = df[df['age'].between(0, MAX_AGE)].copy()
    df['age'] = 70  # force AG70 group

    to_day = lambda col: (df[col] - START_DATE).dt.days
    df['death_day'] = to_day('datumumrti')
    for col in dose_cols_lower:
        df[col + '_day'] = to_day(col)

    dose_day_cols = [col + '_day' for col in dose_cols_lower]
    df['first_dose_day'] = df[dose_day_cols].min(axis=1, skipna=True)
    df['has_any_dose'] = df[dose_day_cols].notna().any(axis=1)

    END_MEASURE = int(df['death_day'].dropna().max())
    df['end_day'] = df['death_day'].fillna(END_MEASURE)

    tv_data = []
    for idx, row in df.iterrows():
        pid = idx
        death_day = row['death_day']
        end_day = row['end_day']
        dose_day = row['first_dose_day']
        has_dose = row['has_any_dose']

        unvax_stop = min(end_day, dose_day + LAG_DAYS) if has_dose else end_day
        tv_data.append({
            'id': pid,
            'start': 0,
            'stop': unvax_stop,
            'event': int(death_day == unvax_stop),
            'vaccinated': 0
        })

        if has_dose and dose_day + LAG_DAYS < end_day:
            tv_data.append({
                'id': pid,
                'start': dose_day + LAG_DAYS,
                'stop': end_day,
                'event': int(death_day == end_day),
                'vaccinated': 1
            })

    tv_df = pd.DataFrame(tv_data)
    tv_df.loc[(tv_df["start"] == tv_df["stop"]) & (tv_df["event"] == 1), "stop"] += 0.5
    tv_df['duration'] = tv_df['stop'] - tv_df['start']

    return df, tv_df, END_MEASURE

def fit_km(tv_df):
    kmf_uvx = KaplanMeierFitter()
    kmf_uvx.fit(tv_df.loc[tv_df['vaccinated'] == 0, 'duration'],
                event_observed=tv_df.loc[tv_df['vaccinated'] == 0, 'event'],
                label="Unvaccinated")

    kmf_vx = KaplanMeierFitter()
    kmf_vx.fit(tv_df.loc[tv_df['vaccinated'] == 1, 'duration'],
               event_observed=tv_df.loc[tv_df['vaccinated'] == 1, 'event'],
               label="Vaccinated")

    return kmf_uvx, kmf_vx

def compute_daily_death_rate_diff(kmf_uvx, kmf_vx, max_day):
    full_index = pd.Index(range(int(max_day)+1))
    surv_uvx = kmf_uvx.survival_function_.reindex(full_index, method='ffill').fillna(1.0)
    surv_vx = kmf_vx.survival_function_.reindex(full_index, method='ffill').fillna(1.0)

    death_rate_uvx = surv_uvx['Unvaccinated'].values[:-1] - surv_uvx['Unvaccinated'].values[1:]
    death_rate_vx = surv_vx['Vaccinated'].values[:-1] - surv_vx['Vaccinated'].values[1:]
    death_rate_uvx = np.append(death_rate_uvx, 0)
    death_rate_vx = np.append(death_rate_vx, 0)

    return full_index.values, death_rate_vx - death_rate_uvx

def compute_daily_dose_counts(df, end_day):
    first_dose_counts = df['first_dose_day'].value_counts().reindex(range(end_day+1), fill_value=0).sort_index()
    all_dose_days = pd.Series(np.concatenate([
        df[f'datum_{i}_day'].dropna().values for i in range(1, 8)
    ])).astype(int)
    all_dose_counts = all_dose_days.value_counts().reindex(range(end_day+1), fill_value=0).sort_index()
    return first_dose_counts, all_dose_counts

# Process simulated data
df_sim, tv_df_sim, end_sim = preprocess_data(SIM_CSV)
kmf_uvx_sim, kmf_vx_sim = fit_km(tv_df_sim)
days_sim, diff_sim = compute_daily_death_rate_diff(kmf_uvx_sim, kmf_vx_sim, end_sim)

# Process real data
df_real, tv_df_real, end_real = preprocess_data(REAL_CSV)
kmf_uvx_real, kmf_vx_real = fit_km(tv_df_real)
days_real, diff_real = compute_daily_death_rate_diff(kmf_uvx_real, kmf_vx_real, end_real)

# Align and compute adjusted diff
max_day = min(end_sim, end_real)
days_common = np.arange(max_day + 1)
diff_sim_aligned = pd.Series(diff_sim, index=days_sim).reindex(days_common, fill_value=0).values
diff_real_aligned = pd.Series(diff_real, index=days_real).reindex(days_common, fill_value=0).values
diff_adjusted = diff_real_aligned - diff_sim_aligned

# Smooth curves
sigma = 3
diff_sim_smooth = gaussian_filter1d(diff_sim_aligned, sigma=sigma)
diff_real_smooth = gaussian_filter1d(diff_real_aligned, sigma=sigma)
diff_adjusted_smooth = gaussian_filter1d(diff_adjusted, sigma=sigma)

# Dose counts
first_dose_real, all_dose_real = compute_daily_dose_counts(df_real, max_day)
# Find vaccination start day (first non-zero day in first dose counts)
vax_start_day = first_dose_real[first_dose_real > 0].index.min()


# Plot
fig = go.Figure()

# Death rate diffs
fig.add_trace(go.Scatter(x=days_common, y=diff_sim_aligned, name='Bias Baseline (Simulated)', line=dict(color='gray')))
fig.add_trace(go.Scatter(x=days_common, y=diff_real_aligned, name='Observed Effect (Real Data)', line=dict(color='blue')))
fig.add_trace(go.Scatter(x=days_common, y=diff_adjusted, name='Adjusted Effect', line=dict(color='green')))
fig.add_trace(go.Scatter(x=days_common, y=diff_sim_smooth, name='Bias Baseline Smoothed', line=dict(width=0.8, color='gray', dash='solid')))
fig.add_trace(go.Scatter(x=days_common, y=diff_real_smooth, name='Observed Effect Smoothed', line=dict(width=0.8, color='blue', dash='solid')))
fig.add_trace(go.Scatter(x=days_common, y=diff_adjusted_smooth, name='Adjusted Effect Smoothed', line=dict(width=0.8, color='green', dash='solid')))

# Smooth dose counts
sigma_dose = 3
first_dose_smooth = gaussian_filter1d(first_dose_real.values, sigma=sigma_dose)
all_dose_smooth = gaussian_filter1d(all_dose_real.values, sigma=sigma_dose)

# Dose count traces (secondary y-axis) - smoothed
fig.add_trace(go.Scatter(x=days_common, y=first_dose_smooth, name='First Doses (smoothed)', yaxis='y2', line=dict(width=0.8, color='orangered', dash='solid')))
fig.add_trace(go.Scatter(x=days_common, y=all_dose_smooth, name='All Doses (smoothed)', yaxis='y2', line=dict(width=0.8, color='orange', dash='solid')))

# Red dotted vertical line at vax start
fig.add_trace(go.Scatter(
    x=[vax_start_day, vax_start_day],
    y=[min(diff_adjusted_smooth), max(diff_adjusted_smooth)],
    mode='lines',
    name='Vaccination Start',
    line=dict(color='red', width=2, dash='dot'),
    hoverinfo='x+name'
))

fig.update_layout(
    title="Daily Death Rate Difference (Vaccinated - Unvaccinated) for Age 70",
    xaxis_title="Days since Jan 1, 2020",
    yaxis=dict(title="Death Rate Difference"),
    yaxis2=dict(title="Dose Counts", overlaying='y', side='right', showgrid=False),
    template='plotly_white',
    hovermode='x unified',
    #legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    legend=dict(orientation="v",yanchor="top",y=1,xanchor="left",x=1.02,bgcolor='rgba(255,255,255,0.8)',bordercolor='lightgray',borderwidth=1)
)
fig.add_annotation(
    text="Dose counts are plotted on the secondary y-axis on the right.<br>Both first doses per day and all doses per day are included.",
    xref="paper", yref="paper", x=0.5, y=-0.2, showarrow=False, align="center"
)
fig.write_html(OUTPUT_HTML)
print(f"Plot saved to {OUTPUT_HTML}")

# Kaplan-Meier plot
fig_surv = go.Figure()
fig_surv.add_trace(go.Scatter(x=kmf_uvx_real.survival_function_.index, y=kmf_uvx_real.survival_function_['Unvaccinated'], name='Unvaccinated', line=dict(color='red')))
fig_surv.add_trace(go.Scatter(x=kmf_vx_real.survival_function_.index, y=kmf_vx_real.survival_function_['Vaccinated'], name='Vaccinated', line=dict(color='blue')))
fig_surv.update_layout(
    title="Kaplan-Meier Survival Curves for Age 70 (Real Data)",
    xaxis_title="Days",
    yaxis_title="Survival Probability",
    template="plotly_white",
    hovermode="x unified",
    yaxis=dict(range=[0, 1]),
)
fig_surv.write_html(OUTPUT_SURV_DIFF_HTML)
print(f"Plot saved to {OUTPUT_SURV_DIFF_HTML}")
