#!/usr/bin/env python3
import os, time
import cdsapi

OUTDIR = "./ESACCI_1992-2020"
os.makedirs(OUTDIR, exist_ok=True)

YEARS   = range(2016, 2023)   # 2016..2022
VERSION = "v2_1_1"
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
        "version": [VERSION]
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

