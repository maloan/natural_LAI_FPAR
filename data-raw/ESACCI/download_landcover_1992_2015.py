#!/usr/bin/env python3
"""Download ESA CCI land-cover archives for 1992-2015 with retry logic."""
import os
import time

import cdsapi

OUTDIR = "./ESACCI_1992-2020"
os.makedirs(OUTDIR, exist_ok=True)

YEARS   = range(1992, 2016)   # 1992..2015
VERSION = "v2_0_7cds"
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

