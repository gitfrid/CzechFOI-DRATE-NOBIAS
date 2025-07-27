import pandas as pd
import numpy as np
import plotly.graph_objs as go


# This script processes simulated or real-world COVID-19 vaccination and death data
# from the Czech dataset to compute and visualize death rates among vaccinated
# (vx) and unvaccinated (uvx) individuals across ages 0–113 and days since Jan 1, 2020.

# Features:
#  - Supports multiple case scenarios (real/simulated deaths and dose schedules)
#  - Classifies population by vaccination status and computes population sizes and deaths
#  - Calculates raw and normalized (per 100k) death rates by group (vx, uvx, total)
#  - Smooths time series using a rolling mean
#  - Computes differences in death rates between vaccinated and unvaccinated
#  - Aggregates and smooths first and all vaccine dose events
#  - Generates a multi-trace interactive Plotly HTML visualization

# Required Inputs:
# - INPUT_CSV`: CSV file with birth year, death date, and up to 7 vaccine dose dates
# - OUTPUT_HTML`: Path to save the resulting interactive visualization



# 6 input CSV files (replace paths if needed)


#INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\FG) case3_sim_deaths_sim_real_doses_with_constraint.csv"
#OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\ZI) vx uvx norm\ZI-FG) case3_sim_deaths_sim_real_doses_with_constraint AG70.html"

INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_HTML = r"C:\CzechFOI-DRATE-NOBIAS\Plot Results\ZI) vx uvx norm\ZI) real data Vesely_106_202403141131_AG70.html"


START_DATE = pd.Timestamp('2020-01-01')
MAX_AGE = 113
REFERENCE_YEAR = 2023

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

# === Simulation Time Frame and Data Structures ===
END_MEASURE = int(df['death_day'].dropna().max())
days = np.arange(0, END_MEASURE + 1)
ages = np.arange(0, MAX_AGE + 1)

results = {
    'day': [],
    'age': [],
    'pop_vx': [],
    'pop_uvx': [],
    'death_vx': [],
    'death_uvx': [],
    'death_total': [],
    'pop_total': [],
}

df_age_groups = [df[df['age'] == age] for age in ages]

# === Main Loop ===
for age, sub in zip(ages, df_age_groups):
    if sub.empty:
        continue

    death_days = sub['death_day'].values
    first_dose_days = sub['first_dose_day'].values
    has_any_dose = sub['has_any_dose'].values

    for day in days:
        alive_mask = np.isnan(death_days) | (death_days > day)
        death_today_mask = (death_days == day)

        is_vaxed = (day >= first_dose_days) & has_any_dose
        is_uvx = ~is_vaxed

        pop_vx = np.sum(alive_mask & is_vaxed)
        pop_uvx = np.sum(alive_mask & is_uvx)
        pop_total = pop_vx + pop_uvx

        death_vx = np.sum(death_today_mask & is_vaxed)
        death_uvx = np.sum(death_today_mask & is_uvx)
        death_total = death_vx + death_uvx

        results['day'].append(day)
        results['age'].append(age)
        results['pop_vx'].append(pop_vx)
        results['pop_uvx'].append(pop_uvx)
        results['pop_total'].append(pop_total)
        results['death_vx'].append(death_vx)
        results['death_uvx'].append(death_uvx)
        results['death_total'].append(death_total)

# === Normalize and Smooth ===
result_df = pd.DataFrame(results)
result_df['deathdiff_uvx_vx'] = result_df['death_uvx'] - result_df['death_vx']
result_df['death_vx_norm'] = (result_df['death_vx'] / result_df['pop_vx'].replace(0, np.nan)) * 100_000
result_df['death_uvx_norm'] = (result_df['death_uvx'] / result_df['pop_uvx'].replace(0, np.nan)) * 100_000
result_df['death_total_norm'] = (result_df['death_total'] / result_df['pop_total'].replace(0, np.nan)) * 100_000
result_df['deathdiff_uvx_vx_norm'] = result_df['death_uvx_norm'] - result_df['death_vx_norm']
# Optional: normalized death difference
result_df.fillna(0, inplace=True)

window_size = 7
result_df['death_vx_norm_smooth'] = result_df.groupby('age')['death_vx_norm'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['death_uvx_norm_smooth'] = result_df.groupby('age')['death_uvx_norm'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['death_total_norm_smooth'] = result_df.groupby('age')['death_total_norm'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['deathdiff_uvx_vx_norm_smooth'] = result_df.groupby('age')['deathdiff_uvx_vx_norm'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['death_vx_smooth'] = result_df.groupby('age')['death_vx'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['death_uvx_smooth'] = result_df.groupby('age')['death_uvx'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['death_total_smooth'] = result_df.groupby('age')['death_total'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)
result_df['deathdiff_uvx_vx_smooth'] = result_df.groupby('age')['deathdiff_uvx_vx'].transform(
    lambda x: x.rolling(window_size, center=True, min_periods=1).mean()
)




# === Dose Counts ===
first_dose_counts_age = {age: pd.Series(0, index=days, dtype=float) for age in ages}
all_dose_counts_age = {age: pd.Series(0, index=days, dtype=float) for age in ages}

for age, sub in zip(ages, df_age_groups):
    if sub.empty:
        continue

    first_counts = sub['first_dose_day'].value_counts().dropna().astype(int)
    s_first = pd.Series(0, index=days, dtype=float)
    s_first.update(first_counts)
    first_dose_counts_age[age] = s_first

    all_dose_days = pd.concat([sub[col + '_day'] for col in dose_date_cols_lower])
    all_counts = all_dose_days.value_counts().dropna().astype(int)
    s_all = pd.Series(0, index=days, dtype=float)
    s_all.update(all_counts)
    all_dose_counts_age[age] = s_all

first_dose_df = pd.DataFrame(first_dose_counts_age)
all_dose_df = pd.DataFrame(all_dose_counts_age)

first_dose_df_smooth = first_dose_df.rolling(window=window_size, center=True, min_periods=1).mean()
all_dose_df_smooth = all_dose_df.rolling(window=window_size, center=True, min_periods=1).mean()

# === Plotly Visualization ===
fig = go.Figure()
colors_vx = 'rgba(0,100,255,0.3)'
colors_uvx = 'rgba(255,0,0,0.3)'
colors_total = 'rgba(0,0,0,0.3)'
colors_diff = 'rgba(0,200,0,0.5)'  # greenish for difference

for age in ages:
    df_age = result_df[result_df['age'] == age]
    if df_age.empty:
        continue
    
    # Norm smooth
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_vx_norm_smooth'],
                             name=f'death_vx_norm_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_vx), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_uvx_norm_smooth'],
                             name=f'death_uvx_norm_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_uvx), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_total_norm_smooth'],
                             name=f'death_total_norm_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_total), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['deathdiff_uvx_vx_norm_smooth'],
                            name=f'deathdiff_uvx_vx_norm_smooth age {age}', yaxis='y1',
                            mode='lines', line=dict(width=1, color=colors_diff), visible='legendonly')) 
    # Norm raw
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_vx_norm'],
                             name=f'death_vx_norm age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_vx.replace('0.3', '0.5')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_uvx_norm'],
                             name=f'death_uvx_norm age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_uvx.replace('0.3', '0.5')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_total_norm'],
                             name=f'death_total_norm age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_total.replace('0.3', '0.5')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['deathdiff_uvx_vx_norm'],
                            name=f'deathdiff_uvx_vx_norm age {age}', yaxis='y1',
                            mode='lines', line=dict(width=1, color=colors_total.replace('0.3', '0.5')), visible='legendonly'))
    # Raw death counts smooth
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_vx_smooth'],
                             name=f'death_vx_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_vx), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_uvx_smooth'],
                             name=f'death_uvx_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_uvx), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_total_smooth'],
                             name=f'death_total_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_total), visible='legendonly'))   
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['deathdiff_uvx_vx_smooth'],
                             name=f'deathdiff_uvx_vx_smooth age {age}', yaxis='y1',
                             mode='lines', line=dict(width=1, color=colors_total), visible='legendonly'))   
    # Raw death counts
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_vx'],
                             name=f'death_vx age {age}', yaxis='y2',
                             mode='lines', line=dict(width=1, color=colors_vx.replace('0.3', '0.15')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_uvx'],
                             name=f'death_uvx age {age}', yaxis='y2',
                             mode='lines', line=dict(width=1, color=colors_uvx.replace('0.3', '0.15')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['death_total'],
                             name=f'death_total age {age}', yaxis='y2',
                             mode='lines', line=dict(width=1, color=colors_total.replace('0.3', '0.15')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['deathdiff_uvx_vx'],
                            name=f'death_vx - death_uvx age {age}', yaxis='y2',
                              mode='lines', line=dict(width=1, color=colors_total.replace('0.3', '0.15')), visible='legendonly'))
    # Population
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['pop_vx'],
                             name=f'pop_vx age {age}', yaxis='y3',
                             mode='lines', line=dict(width=1.5, color=colors_vx.replace('0.3', '0.1')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['pop_uvx'],
                             name=f'pop_uvx age {age}', yaxis='y3',
                             mode='lines', line=dict(width=1.5, color=colors_uvx.replace('0.3', '0.1')), visible='legendonly'))
    fig.add_trace(go.Scatter(x=df_age['day'], y=df_age['pop_total'],
                             name=f'pop_total age {age}', yaxis='y3',
                             mode='lines', line=dict(width=1.5, color=colors_total.replace('0.3', '0.1')), visible='legendonly'))

    # Dose counts
    fig.add_trace(go.Scatter(x=days, y=first_dose_df_smooth[age],
                             name=f'First Dose Count (7-day rolling) age {age}', yaxis='y4',
                             mode='lines', line=dict(width=1.5, color='green'), visible='legendonly'))
    fig.add_trace(go.Scatter(x=days, y=all_dose_df_smooth[age],
                             name=f'All Doses Count (7-day rolling) age {age}', yaxis='y4',
                             mode='lines', line=dict(width=1.5, color='orange'), visible='legendonly'))
        
    

# === Layout with multiple y-axes ===
fig.update_layout(
    title='Vaccinated vs Unvaccinated Deaths, Population, and Doses by Age (Timeline Based on Deaths Only)',
    xaxis=dict(title='Days since 2020-01-01'),
    yaxis=dict(title='Normalized Death/Deathdiff Rate per 100k', side='left', autorange=True),
    yaxis2=dict(title='Raw Deaths/Raw Deatdiff  ', overlaying='y', side='right', position=0.95, autorange=True),
    yaxis3=dict(title='Population', overlaying='y', side='right', position=1.0, autorange=True), #, type='log'
    yaxis4=dict(title='Dose Counts (7-day rolling)', overlaying='y', side='left', position=0.05, autorange=True),
    template='plotly_white',
    height=900,
    showlegend=True
)

fig.write_html(OUTPUT_HTML)
print(f"Plot saved to {OUTPUT_HTML}")
