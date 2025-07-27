import pandas as pd
import numpy as np
import os

# === CONFIGURABLE CONSTANTS ===
INPUT_CSV = r"C:\CzechFOI-DRATE-NOBIAS\Terra\Vesely_106_202403141131_AG70.csv"
OUTPUT_FOLDER = r"C:\CzechFOI-DRATE-NOBIAS\Terra"

START_DATE = pd.Timestamp('2020-01-01')
DOSE_DATE_COLS = [f'Datum_{i}' for i in range(1, 8)]
NEEDED_COLS = ['Rok_narozeni', 'DatumUmrti'] + DOSE_DATE_COLS

RETRIES = 10000
BASE_RNG_SEED = 42

np.random.seed(BASE_RNG_SEED)

# === UTILITIES ===

def to_day_number(date_series):
    date_series = pd.to_datetime(date_series, errors='coerce')
    return (date_series - START_DATE).dt.days

def parse_dates(df):
    for col in DOSE_DATE_COLS + ['DatumUmrti']:
        df[col] = pd.to_datetime(df[col], errors='coerce')
    return df

def estimate_death_rate(df):
    deaths = df["DatumUmrti"]
    death_rate = np.clip(deaths.notna().sum() / len(deaths), 1e-4, 0.999)
    return death_rate

def simulate_deaths(df, end_measure, death_rate):
    df = df.copy()
    df['DatumUmrti'] = pd.NaT
    df['death_day'] = np.nan

    n = len(df)
    will_die = np.random.rand(n) < death_rate
    death_days = np.full(n, np.nan)
    death_days[will_die] = np.random.randint(0, end_measure + 1, size=will_die.sum())

    df.loc[:, "death_day"] = death_days
    df.loc[will_die, "DatumUmrti"] = START_DATE + pd.to_timedelta(death_days[will_die], unit='D')

    return df

# === DOSE ASSIGNMENT MERGED ===

def assign_doses_real_curve_random(df_target, df_source, retries=10000):
    df_target = df_target.copy()
    df_target[DOSE_DATE_COLS] = pd.NaT

    dose_sets = df_source[DOSE_DATE_COLS].dropna(how='all').values.tolist()
    death_day_arr = df_target["death_day"].to_numpy()
    vax_stat_arr = np.zeros(len(death_day_arr), dtype=np.int8)
    rng = np.random.default_rng(BASE_RNG_SEED)

    updates = []
    skip_count = 0

    for dose_dates in dose_sets:
        valid_dates = [d for d in dose_dates if pd.notna(d)]
        if not valid_dates:
            continue

        valid_days = np.array(to_day_number(pd.Series(valid_dates)))
        last_dose_day = valid_days.max()

        eligible_indices = np.where(vax_stat_arr == 0)[0]
        if eligible_indices.size == 0:
            skip_count += 1
            continue

        rng.shuffle(eligible_indices)
        trial_pool = rng.choice(eligible_indices, size=min(retries, len(eligible_indices)), replace=False)

        selected_pos = None
        for trial_pos in trial_pool:
            if np.isnan(death_day_arr[trial_pos]) or death_day_arr[trial_pos] > last_dose_day:
                selected_pos = trial_pos
                break

        if selected_pos is not None:
            updates.append((selected_pos, dose_dates))
            vax_stat_arr[selected_pos] = 1
        else:
            skip_count += 1

    for pos, dose_dates in updates:
        for j, d in enumerate(dose_dates):
            if pd.notna(d):
                df_target.iloc[pos, df_target.columns.get_loc(DOSE_DATE_COLS[j])] = pd.Timestamp(d)

    print(f"Assigned {len(updates)} doses, Skipped {skip_count})")
    return df_target

# === OUTPUT ===

def format_and_save(df, out_path):
    for col in DOSE_DATE_COLS + ['DatumUmrti']:
        df[col] = pd.to_datetime(df[col], errors='coerce').dt.strftime('%Y-%m-%d').fillna('')
    df.to_csv(out_path, index=False)

def save_case(df, filename):
    out_path = os.path.join(OUTPUT_FOLDER, filename)
    format_and_save(df, out_path)
    print(f"Saved: {out_path}")
    return out_path

# === MAIN ===

def run_all_cases():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    print("📥 Loading data...")
    df = pd.read_csv(INPUT_CSV, usecols=NEEDED_COLS, dtype=str)
    df = parse_dates(df)

    max_death_day = to_day_number(df["DatumUmrti"]).max()
    END_MEASURE = max_death_day if not np.isnan(max_death_day) else 1533
    print(f"Measurement window (END_MEASURE): {END_MEASURE} days")

    death_rate = estimate_death_rate(df)
    df_sim_deaths = simulate_deaths(df, end_measure=END_MEASURE, death_rate=death_rate)

    # Case 3: Sim deaths, simulated doses with constraint
    df_case3 = assign_doses_real_curve_random(df_sim_deaths.copy(), df)
    save_case(df_case3, "FG) case3_sim_deaths_sim_real_doses_with_constraint.csv")

    print("✅ All cases processed and saved.")

if __name__ == "__main__":
    run_all_cases()
