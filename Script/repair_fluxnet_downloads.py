"""
UTILITY (not a pipeline step) - Repair broken FLUXNET downloads.

Step 11 skips every site that already has a file in data_raw/, so a download
that was cut off (0-byte or truncated zip) is never retried and the site
silently drops out of the whole analysis. This script
  1. tests every zip in data_raw/ and moves the unreadable ones to
     data_raw/_corrupt/ (moved, not deleted),
  2. downloads those sites again with fluxnet-shuttle, a few at a time,
     checks each new zip, and retries the ones that still fail.
It does NOT rebuild data/fluxnet_daily_all_vars.csv: run 11_fluxnet_download.py
afterwards (it also fetches target sites that have no file at all), then
steps 12-23 for the added sites (step 12 resumes and only extracts new sites).

Usage : python repair_fluxnet_downloads.py            repair everything
        python repair_fluxnet_downloads.py CH-Dav CH-Lae   only these sites
"""
import glob
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path

import pandas as pd
from fluxnet_shuttle import download, listall

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
RAW_DIR = ROOT / "data_raw"
CORRUPT_DIR = RAW_DIR / "_corrupt"
BATCH_SIZE = 4
MAX_ATTEMPTS = 3
SITE_RE = re.compile(r'[A-Z]{2}-[A-Za-z0-9]{3}')


def zip_ok(path):
    try:
        with zipfile.ZipFile(path) as z:
            return z.testzip() is None
    except (zipfile.BadZipFile, OSError):
        return False


def site_of(path):
    m = SITE_RE.search(os.path.basename(path))
    return m.group(0) if m else None


def usable_snapshot():
    """A fresh site catalogue, or the newest saved one that is not empty."""
    try:
        snap = listall(output_dir=str(DATA_DIR))
        if len(pd.read_csv(snap)) > 100:
            return snap
    except Exception as e:
        print(f"  could not refresh the site catalogue ({e}); using a saved one")
    for snap in sorted(glob.glob(str(DATA_DIR / "fluxnet_shuttle_snapshot_*.csv")), reverse=True):
        if len(pd.read_csv(snap)) > 100:
            return snap
    raise SystemExit("No usable site catalogue (snapshot) found.")


only = set(sys.argv[1:])
zips = sorted(glob.glob(str(RAW_DIR / "*.zip")))
print(f"Testing {len(zips)} zip files in '{RAW_DIR}' ...", flush=True)
bad = [z for z in zips if (not only or site_of(z) in only) and not zip_ok(z)]
todo = sorted({site_of(z) for z in bad} | (only - {site_of(z) for z in zips}))
print(f"{len(bad)} unreadable zip(s); {len(todo)} site(s) to download: {' '.join(todo)}", flush=True)
if not todo:
    raise SystemExit("Nothing to repair.")

CORRUPT_DIR.mkdir(exist_ok=True)
for z in bad:
    shutil.move(z, CORRUPT_DIR / os.path.basename(z))
print(f"Moved {len(bad)} file(s) to '{CORRUPT_DIR}'.", flush=True)

snapshot = usable_snapshot()
print(f"Site catalogue: {snapshot}", flush=True)

pending, done = list(todo), []
for attempt in range(1, MAX_ATTEMPTS + 1):
    if not pending:
        break
    print(f"\nAttempt {attempt}: {len(pending)} site(s)", flush=True)
    still = []
    for i in range(0, len(pending), BATCH_SIZE):
        batch = pending[i:i + BATCH_SIZE]
        try:
            download(site_ids=batch, snapshot_file=snapshot, output_dir=str(RAW_DIR))
        except Exception as e:
            print(f"  batch {batch}: {type(e).__name__}: {str(e)[:200]}", flush=True)
        for s in batch:
            files = [z for z in glob.glob(str(RAW_DIR / "*.zip")) if site_of(z) == s]
            if files and all(zip_ok(z) for z in files):
                done.append(s)
                print(f"  OK   {s}  ({sum(os.path.getsize(z) for z in files) / 1e6:.0f} MB)", flush=True)
            else:
                for z in files:                       # incomplete again: set aside so the retry starts clean
                    shutil.move(z, CORRUPT_DIR / (os.path.basename(z) + f".attempt{attempt}"))
                still.append(s)
                print(f"  FAIL {s}", flush=True)
    pending = still

print(f"\nRepaired {len(done)} of {len(todo)} sites.")
if pending:
    print(f"Still failing after {MAX_ATTEMPTS} attempts: {' '.join(pending)}")
print("Next: python Script/11_fluxnet_download.py  (rebuilds the daily file), then steps 12-23.")
