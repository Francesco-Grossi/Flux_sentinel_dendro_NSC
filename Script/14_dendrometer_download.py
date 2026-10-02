"""
PIPELINE STEP 14 - Download the open dendrometer (stem growth) datasets that
cover flux sites. Step 29 turns them into growth per site-year.

    dendrought2018/     Salomon et al. 2022, Zenodo 5711706 (CC-BY 4.0): daily
                        growth and tree water deficit, 77 European sites,
                        2015-2018 (CH-Dav, CZ-BK1, CZ-RAJ, CZ-Stn, FR-Fon)
    harvard_forest/     Harvard Forest Data Archive HF069 (CC0): dendrometer
                        bands on all trees of the EMS tower plots, since 1993
                        (US-Ha1); HF149: Hemlock and LPH towers
    zoebelboden_lter/   LTER Zoebelboden, Zenodo 22940658 (CC-BY 4.0): monthly
                        stem growth of spruce and beech, 1996-2019 (AT-Zoe)
    hyytiala_liu2023/   Liu et al., Zenodo 10037224 (CC-BY 4.0): 30-minute
                        stem radius of two Scots pines, 2015-2019 (FI-Hyy)

Files already present are kept (delete one to refresh it). Datasets that
cannot be fetched automatically (TreeNet, the Dryad oak dataset with US-Ton)
are listed in Output/dendrometer_datasets.md.

Output: data_raw/dendro/<dataset>/...
"""
import re
import zipfile
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data_raw" / "dendro"

ZENODO = {                       # folder -> (record id, files to take; None = all)
    'dendrought2018': (5711706, ['Data - Daily aggregation.xlsx', 'Data - Metadata sites.xlsx']),
    'zoebelboden_lter': (22940658, None),
    'hyytiala_liu2023': (10037224, None),
}
HARVARD_EML = "https://harvardforest1.fas.harvard.edu/data/eml/{}.xml"
HARVARD = {'hf069': ['hf069-09-ems-trees.csv'], 'hf149': None}     # record -> files (None = all csv)


def fetch(url, path):
    if path.exists():
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + '.part')
    with requests.get(url, stream=True, timeout=900) as r:
        r.raise_for_status()
        with open(part, 'wb') as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    part.replace(path)
    print(f"  {path.relative_to(ROOT)} ({path.stat().st_size / 1e6:.1f} MB)")
    return True


n_new = 0
for folder, (record, wanted) in ZENODO.items():
    files = requests.get(f"https://zenodo.org/api/records/{record}", timeout=120).json()['files']
    for f in files:
        if wanted is not None and f['key'] not in wanted:
            continue
        path = OUT / folder / f['key']
        if fetch(f['links']['self'], path):
            n_new += 1
            if path.suffix == '.zip':
                with zipfile.ZipFile(path) as z:
                    z.extractall(path.parent)

for record, wanted in HARVARD.items():
    eml = requests.get(HARVARD_EML.format(record), timeout=120)
    eml.raise_for_status()
    for url in sorted(set(re.findall(r'https?://[^<"\s]+\.csv', eml.text))):
        name = url.rsplit('/', 1)[-1]
        if wanted is None or name in wanted:
            n_new += fetch(url, OUT / 'harvard_forest' / name)

print(f"{n_new} file(s) downloaded; dendrometer data in '{OUT}'.")
