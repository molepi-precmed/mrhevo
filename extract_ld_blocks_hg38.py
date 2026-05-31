#!/usr/bin/env python3
"""
extract_ld_blocks_hg38.py
Read UKBB EUR zarr LD block boundaries (GRCh37) from genoscores via SSH,
liftover each boundary position to hg38, and save the result as JSON.

Output: refpop/ld_blocks_hg38.json
  {"1": [[start_hg38, end_hg38], ...], "2": [...], ..., "22": [...]}
  Each entry is a list of [start_bp, end_bp] pairs in hg38 coordinates.
"""
import json, subprocess, sys
import numpy as np
from pyliftover import LiftOver

ZARR_DIR     = '/opt/datastore/genome/LD_Eur_18mvariants/int8'
GENOSCORES   = 'pmckeigue@genoscores.cphs.mvm.ed.ac.uk'
CHAIN_FILE   = 'refpop/hg19ToHg38.over.chain.gz'
OUTPUT_JSON  = 'refpop/ld_blocks_hg38.json'

# ---- Step 1: fetch .zattrs for all chromosomes from genoscores via SSH ----
print('Fetching block boundaries from genoscores...')
py_cmd = r"""
import json, os
zarr_dir = '/opt/datastore/genome/LD_Eur_18mvariants/int8'
out = {}
for c in range(1, 23):
    path = f'{zarr_dir}/chr_{c}/.zattrs'
    with open(path) as f:
        meta = json.load(f)
    blocks = meta['Estimator properties']['LD blocks']
    out[str(c)] = blocks
print(json.dumps(out))
"""
result = subprocess.run(
    ['ssh', GENOSCORES, f'python3 -c "{py_cmd}"'],
    capture_output=True, text=True, timeout=30
)
if result.returncode != 0:
    sys.exit(f'SSH failed:\n{result.stderr}')

blocks_grch37 = json.loads(result.stdout)
total_raw = sum(len(v) for v in blocks_grch37.values())
print(f'  Retrieved {total_raw} blocks across 22 chromosomes (GRCh37)')

# ---- Step 2: liftover each boundary position GRCh37 → hg38 ----
print(f'Loading chain file {CHAIN_FILE}...')
lo = LiftOver(CHAIN_FILE)

blocks_hg38 = {}
n_lifted = 0
n_merged  = 0

for chrom_str, raw_blocks in sorted(blocks_grch37.items(), key=lambda x: int(x[0])):
    chrom = int(chrom_str)
    lifted = []

    for start37, end37 in raw_blocks:
        # pyliftover is 0-based; our bp positions are 1-based
        r_s = lo.convert_coordinate(f'chr{chrom}', start37 - 1)
        r_e = lo.convert_coordinate(f'chr{chrom}', end37   - 1)

        if r_s and r_e:
            s38 = int(r_s[0][1]) + 1   # back to 1-based
            e38 = int(r_e[0][1]) + 1
            if s38 <= e38:
                lifted.append([s38, e38])
                n_lifted += 1
            else:
                # Positions mapped but inverted (rare, centromeric inversion)
                lifted.append(None)
        else:
            lifted.append(None)   # unmapped

    # Merge None entries into adjacent blocks
    merged = []
    i = 0
    while i < len(lifted):
        if lifted[i] is None:
            n_merged += 1
            # Extend previous block to cover the gap if possible
            if merged:
                # Find next non-None block
                j = i + 1
                while j < len(lifted) and lifted[j] is None:
                    j += 1
                if j < len(lifted):
                    merged[-1][1] = lifted[j][1]
                    i = j + 1
                    continue
            # No previous block: skip
            i += 1
        else:
            merged.append(list(lifted[i]))
            i += 1

    blocks_hg38[chrom_str] = merged

total_out = sum(len(v) for v in blocks_hg38.values())
print(f'  Lifted: {n_lifted} blocks, merged/dropped: {n_merged}')
print(f'  Output: {total_out} blocks across 22 chromosomes (hg38)')

# ---- Step 3: save ----
with open(OUTPUT_JSON, 'w') as f:
    json.dump(blocks_hg38, f, indent=2)
print(f'Saved: {OUTPUT_JSON}')

# Quick check
for c in ['1', '2', '22']:
    first = blocks_hg38[c][0]
    last  = blocks_hg38[c][-1]
    print(f'  chr{c}: {len(blocks_hg38[c])} blocks  '
          f'first=[{first[0]:,}, {first[1]:,}]  last=[{last[0]:,}, {last[1]:,}]')
