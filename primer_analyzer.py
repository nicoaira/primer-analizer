#!/usr/bin/env python3
"""Simple primer analyzer CLI."""

import argparse
import sys

VALID_BASES = set("ACGT")
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Average molecular weights (g/mol) for dNMPs in ssDNA (minus water for internal bonds).
# Using Promega/IDT standard approximation:
#   MW = (nA * 313.21) + (nT * 304.2) + (nC * 289.18) + (nG * 329.21) - 61.96
MW_A, MW_T, MW_C, MW_G, MW_ADJUST = 313.21, 304.2, 289.18, 329.21, -61.96

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998).
# Values are (dH kcal/mol, dS cal/(mol*K)).
NN_PARAMS = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)   # initiation with terminal G/C
INIT_AT = (2.3, 4.1)    # initiation with terminal A/T
R = 1.987  # cal/(mol*K)


def clean(seq: str) -> str:
    return seq.strip().upper().replace(" ", "").replace("\n", "")


def validate(seq: str) -> None:
    bad = set(seq) - VALID_BASES
    if bad:
        raise ValueError(f"Invalid base(s) in sequence: {', '.join(sorted(bad))}")
    if not seq:
        raise ValueError("Sequence is empty.")


def gc_content(seq: str) -> float:
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def molecular_weight(seq: str) -> float:
    return (
        seq.count("A") * MW_A
        + seq.count("T") * MW_T
        + seq.count("C") * MW_C
        + seq.count("G") * MW_G
        + MW_ADJUST
    )


def tm_wallace(seq: str) -> float:
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc


def tm_basic(seq: str) -> float:
    n = len(seq)
    gc = seq.count("G") + seq.count("C")
    return 64.9 + 41 * (gc - 16.4) / n


def tm_nearest_neighbor(seq: str, primer_nM: float, salt_mM: float) -> float:
    import math
    dH = 0.0
    dS = 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS[seq[i:i + 2]]
        dH += h
        dS += s
    for end in (seq[0], seq[-1]):
        h, s = INIT_GC if end in "GC" else INIT_AT
        dH += h
        dS += s
    ct = primer_nM * 1e-9
    tm_kelvin = (dH * 1000) / (dS + R * math.log(ct / 4))
    return tm_kelvin - 273.15 + 16.6 * math.log10(salt_mM / 1000)


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def gc_clamp(seq: str) -> int:
    """Count G/C in the last 5 bases (ideal: 1-2)."""
    return sum(1 for b in seq[-5:] if b in "GC")


def analyze(seq: str, primer_nM: float, salt_mM: float) -> dict:
    n = len(seq)
    stats = {
        "Sequence (5'->3')": seq,
        "Reverse complement": reverse_complement(seq),
        "Length (nt)": n,
        "GC content (%)": round(gc_content(seq), 2),
        "Molecular weight (g/mol)": round(molecular_weight(seq), 2),
        "GC clamp (last 5 nt)": gc_clamp(seq),
        "Tm Wallace (°C)": round(tm_wallace(seq), 2),
        "Tm basic (°C)": round(tm_basic(seq), 2),
        "Tm nearest-neighbor (°C)": round(
            tm_nearest_neighbor(seq, primer_nM, salt_mM), 2
        ),
    }
    if n < 14:
        stats["Recommended Tm"] = f"{stats['Tm Wallace (°C)']} °C (Wallace, short primer)"
    else:
        stats["Recommended Tm"] = (
            f"{stats['Tm nearest-neighbor (°C)']} °C (nearest-neighbor)"
        )
    return stats


def format_report(stats: dict) -> str:
    width = max(len(k) for k in stats)
    return "\n".join(f"{k.ljust(width)} : {v}" for k, v in stats.items())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Analyze a DNA primer sequence and report basic stats."
    )
    parser.add_argument(
        "sequence",
        nargs="?",
        help="Primer sequence (A/C/G/T). If omitted, read from stdin.",
    )
    parser.add_argument(
        "--primer-conc",
        type=float,
        default=250.0,
        help="Primer concentration in nM for nearest-neighbor Tm (default: 250).",
    )
    parser.add_argument(
        "--salt",
        type=float,
        default=50.0,
        help="Monovalent salt [Na+] in mM (default: 50).",
    )
    args = parser.parse_args(argv)

    raw = args.sequence if args.sequence is not None else sys.stdin.read()
    seq = clean(raw)
    try:
        validate(seq)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    stats = analyze(seq, args.primer_conc, args.salt)
    print(format_report(stats))
    return 0


if __name__ == "__main__":
    sys.exit(main())
