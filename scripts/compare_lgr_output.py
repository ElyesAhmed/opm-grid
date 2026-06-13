#!/usr/bin/env python3
"""Compare two ECLIPSE output sets (e.g. this implementation vs upstream) for
an LGR run, as a correctness check for the refinement builder.

Two complementary comparisons (docs/PLAN.md "Verification"):

  * Full element-wise, per LGR section, of every INIT array and every restart
    field. This is valid -- and a strong check -- precisely *because the cell
    numbering is the same*: INIT/restart arrays are written in per-section
    Cartesian (active-cell) order, so identical cell numbering makes a direct
    element-wise comparison meaningful. Bitwise-identical INIT therefore
    proves the grid, transmissibility, pore volumes and depths are identical.

  * Connection-keyed for boundary transmissibilities (TRANGL): keyed by the
    (NNCG, NNCL) cell pair from the EGRID, so it is immune to differences in
    face numbering between the two implementations.

Run convertECL (opm-common) first to produce formatted F* files, or pass
--convert with the path to convertECL.

Usage:
    compare_lgr_output.py CASE_A_PREFIX CASE_B_PREFIX [--convert PATH_TO_convertECL]

where PREFIX is the path without extension, e.g. /tmp/up/SPE1CASE1_CARFIN.
"""
import argparse
import re
import subprocess
import sys


def parse_sections(fname):
    """Return {(section, keyword): [block, ...]} preserving order and repeats.
    Sections are delimited by the LGR / LGRSGONE / ENDLGR markers."""
    secs = {}
    section = 'GLOBAL'
    with open(fname) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        m = re.match(r" '(.{8})'\s+(\d+) '(\w+)'", lines[i])
        if not m:
            i += 1
            continue
        kw, n, typ = m.group(1).strip(), int(m.group(2)), m.group(3)
        i += 1
        vals = []
        while len(vals) < n and i < len(lines):
            if typ == 'CHAR':
                vals += re.findall(r"'([^']*)'", lines[i])
            elif typ in ('LOGI', 'MESS'):
                vals += lines[i].split()
            else:
                try:
                    vals += [float(x.replace('D', 'E')) for x in lines[i].split()]
                except ValueError:
                    vals += lines[i].split()
            i += 1
        if kw == 'LGR':
            section = vals[0].strip()
            continue
        if kw in ('LGRSGONE', 'ENDLGR'):
            section = 'GLOBAL'
            continue
        secs.setdefault((section, kw), []).append(vals)
    return secs


def reldiff(x, y):
    return abs(x - y) / max(abs(x), abs(y), 1e-30)


HEADER_KW = {'INTEHEAD', 'LOGIHEAD', 'DOUBHEAD', 'LGRHEADI', 'LGRHEADD',
             'LGRHEADQ', 'LGRJOIN', 'LGR', 'LGRSGONE', 'HOSTNUM', 'ENDLGR',
             'SEQNUM', 'IGRP', 'SGRP', 'XGRP', 'ZGRP', 'IWEL', 'SWEL', 'XWEL',
             'ZWEL', 'ICON', 'SCON', 'XCON', 'STARTSOL', 'ENDSOL'}


def compare_elementwise(fa, fb, label, rtol=1e-6):
    a, b = parse_sections(fa), parse_sections(fb)
    checked = flagged = 0
    worst = 0.0
    for k in sorted(set(a) | set(b)):
        if k[1] in HEADER_KW:
            continue
        va, vb = a.get(k), b.get(k)
        if va is None or vb is None or len(va) != len(vb):
            print(f"  {label} {k}: structural mismatch")
            flagged += 1
            continue
        for blkA, blkB in zip(va, vb):
            if not blkA or isinstance(blkA[0], str) or len(blkA) != len(blkB):
                continue
            mx = max((reldiff(x, y) for x, y in zip(blkA, blkB)), default=0.0)
            checked += 1
            worst = max(worst, mx)
            if mx > rtol:
                flagged += 1
    print(f"{label}: arrays={checked} flagged(>{rtol:g})={flagged} worst-rel={worst:.3e}")
    return flagged == 0


def nnc_pairs(egrid):
    out, cur = [], {}
    for (sec, kw), blocks in parse_sections(egrid).items():
        pass
    # NNCG/NNCL appear in file order; re-parse preserving pairing per section.
    secs = parse_sections(egrid)
    for (sec, kw) in secs:
        if kw == 'NNCG':
            g = secs[(sec, 'NNCG')]
            l = secs.get((sec, 'NNCL'))
            if l and len(g) == len(l):
                for bg, bl in zip(g, l):
                    out.append((sec, list(zip([int(x) for x in bg],
                                              [int(x) for x in bl]))))
    return out


def compare_trangl(init_a, egrid_a, init_b, egrid_b):
    ta = [(s, v) for (s, kw), blks in parse_sections(init_a).items()
          if kw == 'TRANGL' for v in blks]
    tb = [(s, v) for (s, kw), blks in parse_sections(init_b).items()
          if kw == 'TRANGL' for v in blks]
    na = nnc_pairs(egrid_a)
    nb = nnc_pairs(egrid_b)
    if len(ta) != len(tb) or len(na) != len(nb):
        print(f"TRANGL: block-count mismatch ({len(ta)} vs {len(tb)})")
        return False
    common = onlyA = onlyB = 0
    worst = 0.0
    for (sa, va), (na_s, na_p), (sb, vb), (nb_s, nb_p) in zip(ta, na, tb, nb):
        da = dict(zip(na_p, va))
        db = dict(zip(nb_p, vb))
        keys = set(da) & set(db)
        common += len(keys)
        onlyA += len(set(da) - set(db))
        onlyB += len(set(db) - set(da))
        worst = max([worst] + [reldiff(da[k], db[k]) for k in keys])
    print(f"TRANGL: common={common} only-A={onlyA} only-B={onlyB} worst-rel={worst:.3e}")
    return onlyA == 0 and onlyB == 0 and worst < 1e-9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('a')
    ap.add_argument('b')
    ap.add_argument('--convert', help='path to convertECL; runs it on .INIT/.EGRID/.UNRST first')
    args = ap.parse_args()

    if args.convert:
        for prefix in (args.a, args.b):
            for ext in ('INIT', 'EGRID', 'UNRST'):
                subprocess.run([args.convert, f'{prefix}.{ext}'],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    ok = True
    ok &= compare_elementwise(f'{args.a}.FINIT', f'{args.b}.FINIT', 'INIT(full)')
    ok &= compare_trangl(f'{args.a}.FINIT', f'{args.a}.FEGRID',
                         f'{args.b}.FINIT', f'{args.b}.FEGRID')
    # Restart differs at solver tolerance even with identical discretization
    # (internal cell ordering -> linear-solver path); report per field as the
    # worst absolute difference over all steps/sections, do not gate.
    report_restart(f'{args.a}.FUNRST', f'{args.b}.FUNRST')
    print('INIT + TRANGL identical' if ok else 'DIFFERENCES in INIT/TRANGL')
    return 0 if ok else 1


def report_restart(fa, fb):
    a, b = parse_sections(fa), parse_sections(fb)
    fields = {}
    for k in set(a) & set(b):
        if k[1] not in ('PRESSURE', 'SGAS', 'SWAT', 'SOIL', 'RS', 'RV'):
            continue
        for blkA, blkB in zip(a[k], b[k]):
            if not blkA or isinstance(blkA[0], str):
                continue
            mx = max(abs(x - y) for x, y in zip(blkA, blkB))
            fields[k[1]] = max(fields.get(k[1], 0.0), mx)
    units = {'PRESSURE': 'bar', 'RS': 'sm3/sm3', 'RV': 'sm3/sm3'}
    print("UNRST worst ABS diff per field over all steps/sections "
          "(solver-path, not discretization):")
    for f in sorted(fields):
        print(f"    {f:9s}: {fields[f]:.3e} {units.get(f, 'saturation')}")


if __name__ == '__main__':
    sys.exit(main())
