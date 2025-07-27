import pandas as pd
import os

# === CONFIGURATION ===
INPUT_CSV = r"C:\CzechFOI-DRATE_NOBIAS\TERRA\Vesely_106_202403141131.csv"
OUTPUT_CSV = r"C:\CzechFOI-DRATE_NOBIAS\TERRA\Vesely_106_202403141131_AG70.csv"
REFERENCE_YEAR = 2023
DOSE_DATE_COLS = [f'Datum_{i}' for i in range(1, 8)]
NEEDED_COLS = ['Rok_narozeni', 'DatumUmrti'] + DOSE_DATE_COLS

# === FUNCTIONS ===
def parse_dates(df):
    for col in DOSE_DATE_COLS + ['DatumUmrti']:
        df[col] = pd.to_datetime(df[col], errors='coerce')
    return df

def calculate_age(df):
    df["Age"] = REFERENCE_YEAR - df["Rok_narozeni"].astype(int)
    return df

def format_dates_for_csv(df):
    for col in DOSE_DATE_COLS + ['DatumUmrti']:
        df[col] = df[col].dt.strftime('%Y-%m-%d').fillna('')
    return df

# === MAIN ===
def filter_and_save_age_70():
    print("📥 Loading input CSV...")
    df = pd.read_csv(INPUT_CSV, usecols=NEEDED_COLS, dtype=str)

    print("📆 Parsing dates and calculating age...")
    df = parse_dates(df)
    df = calculate_age(df)

    print("🔎 Filtering to Age == 70...")
    df_ag70 = df[df["Age"] == 70].copy()

    print(f"💾 Saving {len(df_ag70)} rows to output...")
    df_ag70 = format_dates_for_csv(df_ag70)
    df_ag70.to_csv(OUTPUT_CSV, index=False)

    print(f"✅ Done. Saved to {OUTPUT_CSV}")

# === RUN ===
if __name__ == "__main__":
    filter_and_save_age_70()
