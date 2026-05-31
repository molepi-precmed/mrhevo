#!/usr/bin/env python3
"""
run_zarr_loci.py -- compute alpha and/or gamma coefficients for a set of loci
using the UKBB EUR LD zarr store on genoscores.

Reads JSON from stdin, writes JSON to stdout, progress to stderr.

Input JSON schema:
{
  "zarr_dir":     "/opt/datastore/genome/LD_Eur_18mvariants/int8",
  "min_eig_frac": 0.01,
  "N_gamma":      12461,   # omit to skip gamma
  "p_case":       0.398,   # omit to skip gamma
  "loci": [
    {
      "locus_id": "chr1_locus1",
      "chr": 1,
      "snps": [
        { "rsid": "rs1", "Allele1": "A", "Allele2": "G",
          "Effect": 0.01, "Freq1": 0.3,
          "StdErr": 0.001, "TotalSampleSize": 33902,
          "pos": 1234567,
          "z_gamma": 1.23   # optional; include for gamma computation
        }
      ]
    }
  ]
}

Output JSON schema:
{
  "loci": [
    { "locus_id": "chr1_locus1", "chr": 1,
      "locus_start": 1234567, "locus_end": 1234567,
      "n_snps": 5, "eff_rank": 3,
      "alpha_hat": 0.05, "se_alpha_hat": 0.01, "sd_Z": 0.22,
      "gamma_hat": 0.02, "se_gamma_hat": 0.005   # if N_gamma provided
    }
  ]
}
"""
import sys
import json
import numpy as np
import magenpy as mgp
import os

_ld_cache = {}


def _open_chr(zarr_dir, chrom):
    key = (zarr_dir, chrom)
    if key not in _ld_cache:
        path = os.path.join(zarr_dir, f'chr_{chrom}')
        ldm  = mgp.LDMatrix.from_path(path)
        st   = ldm.to_snp_table(col_subset=['SNP', 'POS', 'A1', 'A2', 'MAF'])
        idx  = {rsid: i for i, rsid in enumerate(st['SNP'].values)}
        _ld_cache[key] = (ldm, st, idx)
    return _ld_cache[key]


def _get_ld_submatrix(zarr_dir, chrom, rsids):
    """Return (C_g, ld_snps) or (None, None). ld_snps is a list of dicts."""
    try:
        chrom = int(chrom)
    except (ValueError, TypeError):
        return None, None
    ldm, st, idx = _open_chr(zarr_dir, chrom)
    our_indices = np.array(
        sorted(idx[r] for r in rsids if r in idx), dtype=np.int64
    )
    if len(our_indices) == 0:
        return None, None
    start_row = int(our_indices[0])
    end_row   = int(our_indices[-1]) + 1
    csr = ldm.load_data(
        start_row=start_row, end_row=end_row,
        return_symmetric=True, return_as_csr=True,
        dtype=np.float64
    )
    rel = our_indices - start_row
    C_g = np.asarray(csr[np.ix_(rel, rel)].todense(), dtype=np.float64)
    C_g[np.isnan(C_g)] = 0.0
    matched = st.iloc[our_indices].reset_index(drop=True)
    ld_snps = [
        {'rsid': row['SNP'], 'ld_A1': row['A1'], 'ld_A2': row['A2']}
        for _, row in matched.iterrows()
    ]
    return C_g, ld_snps


def _trunc_pseudoinv(S, min_eig_frac):
    """Return (V_k, d_inv, eff_rank) factored pseudoinverse of symmetric S."""
    eigvals, eigvecs = np.linalg.eigh(S)
    eigvals = eigvals[::-1];  eigvecs = eigvecs[:, ::-1]   # descending
    tol  = eigvals[0] * min_eig_frac
    keep = eigvals > tol
    return eigvecs[:, keep], 1.0 / eigvals[keep], int(keep.sum())


def _align_snps(snps, C_g, ld_snps):
    """
    Merge per-locus SNP data to LD matrix order, aligning Allele1 to LD A1.
    Returns (snps_aln, Cg_sub, keep_idx) or (None, None, None).

    snps_aln: list of dicts with keys Effect, Freq1, StdErr (if present),
              TotalSampleSize (if present), z_gamma (if present), pos.
    """
    snp_map = {s['rsid']: s for s in snps}
    merged   = []
    keep_idx = []
    for i, ld in enumerate(ld_snps):
        rsid = ld['rsid']
        if rsid not in snp_map:
            continue
        s    = snp_map[rsid]
        a1   = s['Allele1'].upper()
        la1  = ld['ld_A1'].upper()
        la2  = ld['ld_A2'].upper()
        if a1 == la1:
            flip = False
        elif a1 == la2:
            flip = True
        else:
            continue   # allele mismatch

        row = dict(s)
        row['Effect'] = -s['Effect'] if flip else s['Effect']
        row['Freq1']  = 1.0 - s['Freq1'] if flip else s['Freq1']
        if 'z_gamma' in s:
            row['z_gamma'] = -s['z_gamma'] if flip else s['z_gamma']
        merged.append(row)
        keep_idx.append(i)

    if not merged:
        return None, None, None

    Cg_sub = C_g[np.ix_(keep_idx, keep_idx)]
    return merged, Cg_sub, keep_idx


def _compute_alpha(snps_aln, Cg, min_eig_frac):
    """Compute alpha coefficients from aligned SNP data."""
    k = len(snps_aln)
    p   = np.array([s['Freq1'] for s in snps_aln])
    eff = np.array([s['Effect'] for s in snps_aln])
    se  = np.array([s['StdErr'] for s in snps_aln])
    N   = np.array([s['TotalSampleSize'] for s in snps_aln])

    sigma_j      = np.sqrt(2.0 * p * (1.0 - p))
    alpha_u_star = eff * sigma_j

    if k == 1:
        alpha_m  = alpha_u_star.copy()
        eff_rank = 1
    else:
        V_k, d_inv, eff_rank = _trunc_pseudoinv(Cg, min_eig_frac)
        alpha_m = V_k @ (d_inv * (V_k.T @ alpha_u_star))

    beta_m    = alpha_m / sigma_j
    alpha_hat = float(np.sqrt(np.sum(beta_m ** 2)))
    if alpha_hat == 0.0:
        return None

    var_S = float(alpha_m @ Cg @ alpha_m)
    if var_S <= 0.0:
        return None
    var_Z = var_S / alpha_hat ** 2

    sigma_X_sq = float(np.median(2.0 * p * (1.0 - p) * (N * se ** 2 + eff ** 2)))
    N_alpha    = float(np.median(N))
    denom      = sigma_X_sq - var_S
    se_alpha   = (1.0 / np.sqrt(N_alpha * var_Z / denom)) if denom > 0.0 else None

    positions = [s['pos'] for s in snps_aln]
    return {
        'n_snps':       k,
        'eff_rank':     eff_rank,
        'locus_start':  int(min(positions)),
        'locus_end':    int(max(positions)),
        'alpha_hat':    alpha_hat,
        'se_alpha_hat': float(se_alpha) if se_alpha is not None else None,
        'sd_Z':         float(np.sqrt(var_Z)),
        # carry through for gamma computation
        '_alpha_m':     alpha_m,
        '_sigma_j':     sigma_j,
        '_var_S':       var_S,
        '_Cg':          Cg,
        '_eff_rank':    eff_rank,
    }


def _compute_gamma(snps_aln, alpha_res, N_gamma, p_case):
    """Compute gamma given aligned SNP data and already-computed alpha internals."""
    gz = np.array([s['z_gamma'] for s in snps_aln])
    gamma_u_star = gz / np.sqrt(N_gamma * p_case * (1.0 - p_case))

    alpha_m   = alpha_res['_alpha_m']
    alpha_hat = alpha_res['alpha_hat']
    var_S     = alpha_res['_var_S']
    var_Z     = var_S / alpha_hat ** 2

    gamma_hat = alpha_hat * float(np.sum(gamma_u_star * alpha_m)) / var_S
    se_gamma  = 1.0 / np.sqrt(N_gamma * var_Z * p_case * (1.0 - p_case))
    return {'gamma_hat': gamma_hat, 'se_gamma_hat': float(se_gamma)}


def process_loci(job):
    zarr_dir     = job['zarr_dir']
    min_eig_frac = job.get('min_eig_frac', 0.01)
    N_gamma      = job.get('N_gamma')
    p_case       = job.get('p_case')
    do_gamma     = (N_gamma is not None) and (p_case is not None)

    results = []
    for locus in job['loci']:
        lid   = locus['locus_id']
        chrom = locus['chr']
        snps  = locus['snps']
        rsids = [s['rsid'] for s in snps]

        C_g, ld_snps = _get_ld_submatrix(zarr_dir, chrom, rsids)
        if C_g is None:
            sys.stderr.write(f'  {lid}: no LD match\n')
            continue

        snps_aln, Cg_sub, _ = _align_snps(snps, C_g, ld_snps)
        if snps_aln is None:
            sys.stderr.write(f'  {lid}: allele alignment failed\n')
            continue

        alpha_res = _compute_alpha(snps_aln, Cg_sub, min_eig_frac)
        if alpha_res is None:
            sys.stderr.write(f'  {lid}: alpha computation failed\n')
            continue

        n_snps   = alpha_res['n_snps']
        eff_rank = alpha_res['eff_rank']
        if eff_rank < n_snps // 10:
            sys.stderr.write(
                f'  {lid}: n_snps={n_snps}, eff_rank={eff_rank} (high LD, truncated)\n'
            )

        r = {
            'locus_id':    lid,
            'chr':         chrom,
            'locus_start': alpha_res['locus_start'],
            'locus_end':   alpha_res['locus_end'],
            'n_snps':      n_snps,
            'eff_rank':    eff_rank,
            'alpha_hat':   alpha_res['alpha_hat'],
            'se_alpha_hat':alpha_res['se_alpha_hat'],
            'sd_Z':        alpha_res['sd_Z'],
        }

        if do_gamma:
            has_zgamma = all('z_gamma' in s for s in snps_aln)
            if has_zgamma:
                gamma_res = _compute_gamma(snps_aln, alpha_res, N_gamma, p_case)
                r.update(gamma_res)
            else:
                sys.stderr.write(f'  {lid}: z_gamma missing, skipping gamma\n')

        results.append(r)
        sys.stderr.write(
            f'  {lid}: n_snps={n_snps}, eff_rank={eff_rank}, '
            f'alpha_hat={alpha_res["alpha_hat"]:.4f}\n'
        )

    # Strip internal keys before returning
    _internal = {'_alpha_m', '_sigma_j', '_var_S', '_Cg', '_eff_rank'}
    for r in results:
        for k in _internal:
            r.pop(k, None)

    return results


def main():
    job     = json.load(sys.stdin)
    results = process_loci(job)
    json.dump({'loci': results}, sys.stdout, indent=2)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
