# Primer Analyzer

A small command-line tool that analyzes a DNA primer sequence and reports basic
stats useful for PCR design: length, GC content, molecular weight, GC clamp,
reverse complement, and three melting-temperature (Tm) estimates.

## What it reports

| Stat | Notes |
|------|-------|
| **Sequence (5'→3')** | The cleaned input (uppercased, whitespace stripped). |
| **Reverse complement** | 5'→3' reverse complement of the input. |
| **Length (nt)** | Number of bases. |
| **GC content (%)** | Fraction of G and C bases. Typical primers: 40–60%. |
| **Molecular weight (g/mol)** | Approximate ssDNA MW. |
| **GC clamp (last 5 nt)** | Count of G/C in the 3' end. Ideal: 1–2. |
| **Tm Wallace (°C)** | `2·(A+T) + 4·(G+C)`. Best for primers ≤14 nt. |
| **Tm basic (°C)** | `64.9 + 41·(G+C−16.4)/N`. Rough estimate for longer primers. |
| **Tm nearest-neighbor (°C)** | SantaLucia 1998 thermodynamics with salt correction. Most accurate. |
| **Recommended Tm** | Wallace for short primers, nearest-neighbor otherwise. |

## Requirements

- Python 3.8 or newer. No third-party dependencies.

## Installation

Clone or copy the file into a directory of your choice:

```bash
git clone <your-repo-url> primer-analizer
cd primer-analizer
```

Optionally make the script directly executable and put it on your `PATH`:

```bash
chmod +x primer_analyzer.py
ln -s "$(pwd)/primer_analyzer.py" ~/.local/bin/primer-analyzer
```

## Usage

Pass the sequence as an argument:

```bash
python3 primer_analyzer.py ATGCGTACGTTAGCTAGCTA
```

Or pipe it from stdin (useful for FASTA-like snippets):

```bash
echo "ATGC GTAC GTTA GCTA GCTA" | python3 primer_analyzer.py
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--primer-conc` | `250` | Primer concentration in nM (used by nearest-neighbor Tm). |
| `--salt` | `50` | Monovalent cation `[Na+]` in mM. |

Example with custom conditions:

```bash
python3 primer_analyzer.py ATGCGTACGTTAGCTAGCTA --primer-conc 500 --salt 75
```

### Example output

```
Sequence (5'->3')        : ATGCGTACGTTAGCTAGCTA
Reverse complement       : TAGCTAGCTAACGTACGCAT
Length (nt)              : 20
GC content (%)           : 45.0
Molecular weight (g/mol) : 6132.06
GC clamp (last 5 nt)     : 2
Tm Wallace (°C)          : 58
Tm basic (°C)            : 49.73
Tm nearest-neighbor (°C) : 46.18
Recommended Tm           : 46.18 °C (nearest-neighbor)
```

## Notes & limitations

- Only accepts unambiguous DNA bases (`A`, `C`, `G`, `T`). Whitespace is ignored;
  lowercase is accepted. Any other character (including `N` or IUPAC codes)
  raises an error.
- Tm values are estimates. Primer design for critical experiments should also
  consider secondary structure, primer-dimers, and specificity — which this
  tool does not check.
- The nearest-neighbor model assumes non-self-complementary duplexes and primer
  in excess over template (`Ct/4` approximation).

## References

- SantaLucia J. (1998). *A unified view of polymer, dumbbell, and oligonucleotide
  DNA nearest-neighbor thermodynamics.* PNAS 95(4):1460–5.
- Wallace R.B. et al. (1979). *Hybridization of synthetic oligodeoxyribonucleotides
  to φX 174 DNA.* Nucleic Acids Research 6(11):3543–57.
