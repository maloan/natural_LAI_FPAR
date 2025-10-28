#!/usr/bin/env python3
import os, time
import cdsapi

OUTDIR = "/home/akurth@giub.local/Documents/SNU/SNU_LAI/data/C3S_LC_2016_2022"
os.makedirs(OUTDIR, exist_ok=True)

YEARS   = range(2016, 2023)   # 2016..2022
VERSION = "v2_1_1"            # <-- underscores!
DATASET = "satellite-land-cover"
MAX_RETRY = 5

c = cdsapi.Client()

for year in YEARS:
    target = os.path.join(OUTDIR, f"C3S_LC_{year}_{VERSION}.zip")
    if os.path.exists(target):
        print(f"[skip] {year}: {target} already exists"); continue

    req = {
        "variable": "all",
        "year": [str(year)],
        "version": [VERSION],
        # "area": [N, W, S, E],   # optional subsetting if you want
        # "format": "zip"         # optional; default is zip for this dataset
    }

    for attempt in range(1, MAX_RETRY + 1):
        try:
            print(f"[{year}] downloading… (attempt {attempt})")
            c.retrieve(DATASET, req, target)
            print(f"[{year}] done")
            break
        except Exception as e:
            if attempt == MAX_RETRY:
                print(f"[{year}] FAILED: {e}")
            else:
                wait = 30 * attempt
                print(f"[{year}] error: {e} — retrying in {wait}s")
                time.sleep(wait)

