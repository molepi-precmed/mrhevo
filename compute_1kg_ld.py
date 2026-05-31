#!/usr/bin/env python3
"""
compute_1kg_ld.py  --  runs on diabepi
Compute a 1000G EUR zarr LD matrix using the same block boundaries as the
UKBB EUR zarr LD matrix (boundaries lifted from GRCh37 to hg38).

Strategy:
  - Load all chromosomes from the combined BED file in one GWADataLoader.
  - For each chromosome, call the GenotypeMatrix.compute_ld() method directly
    with the per-chromosome ld_blocks array.
  - Parallelise over chromosomes using joblib (22 processes).

Usage:
    python3 compute_1kg_ld.py [--chrom N]   # --chrom N for single-chromosome test
    python3 compute_1kg_ld.py               # full run, 22 parallel processes
"""
import argparse, json, os, sys, time, tempfile
import numpy as np
from joblib import Parallel, delayed
import magenpy as mgp

BED_FILE  = os.path.expanduser('~/kg.2020.hg38.eur.bed')
BLOCKS_F  = os.path.expanduser('~/ld_blocks_hg38.json')
OUT_DIR   = os.path.expanduser('~/ld_1kg_eur')
N_JOBS    = 22
# MAF filter (>= 1%) already applied when splitting BED by chromosome with plink

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=int, default=None)
parser.add_argument('--n-jobs', type=int, default=N_JOBS)
args = parser.parse_args()

with open(BLOCKS_F) as f:
    blocks_hg38 = json.load(f)

os.makedirs(OUT_DIR, exist_ok=True)


def compute_one_chr(chrom):
    """Load a single chromosome from the combined BED and compute its LD matrix."""
    t0 = time.time()
    os.environ.update({'OMP_NUM_THREADS': '2', 'MKL_NUM_THREADS': '2',
                       'OPENBLAS_NUM_THREADS': '2'})
    chrom_str = str(chrom)
    if chrom_str not in blocks_hg38:
        return chrom, f'skip: no blocks for chr{chrom}'

    # int32 required by magenpy's map_variants_to_genomic_blocks (pd.merge_asof)
    ld_blocks_arr = np.array(blocks_hg38[chrom_str], dtype=np.int32)

    try:
        # Load per-chromosome BED (pre-split by plink with MAF filter)
        chr_bed = os.path.expanduser(f'~/kg_chr{chrom}.bed')
        if not os.path.exists(chr_bed):
            return chrom, f'skip: {chr_bed} not found'

        loader = mgp.GWADataLoader(bed_files=chr_bed, threads=2)

        if chrom not in loader.genotype:
            return chrom, f'skip: chrom {chrom} not in loader.genotype'

        geno = loader.genotype[chrom]
        n_snps = geno.n_snps

        # Call compute_ld on the per-chromosome GenotypeMatrix directly
        geno.compute_ld(
            estimator='block',
            output_dir=OUT_DIR,
            dtype='int8',
            compressor_name='zstd',
            compression_level=7,
            ld_blocks=ld_blocks_arr
        )

        elapsed = time.time() - t0
        return chrom, f'ok  n_snps={n_snps:,}  n_blocks={len(ld_blocks_arr)}  {elapsed:.0f}s'

    except Exception as e:
        import traceback
        return chrom, f'ERROR: {e}\n{traceback.format_exc()}'


chroms = [args.chrom] if args.chrom else list(range(1, 23))

if len(chroms) == 1:
    c, status = compute_one_chr(chroms[0])
    print(f'chr{c}: {status}')
else:
    print(f'Processing {len(chroms)} chromosomes with up to {args.n_jobs} parallel processes...')
    print(f'Output: {OUT_DIR}')
    t_start = time.time()

    results = Parallel(n_jobs=args.n_jobs, prefer='processes', verbose=5)(
        delayed(compute_one_chr)(c) for c in chroms
    )

    total = time.time() - t_start
    print()
    print('=== Results ===')
    for c, s in sorted(results):
        print(f'  chr{c:2d}: {s}')
    n_ok  = sum(1 for _, s in results if s.startswith('ok'))
    n_err = sum(1 for _, s in results if 'ERROR' in s)
    print(f'\n{n_ok}/22 succeeded, {n_err} failed  ({total:.0f}s total)')
