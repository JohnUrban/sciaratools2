#!/usr/bin/env python3
"""
Convert an IUPAC DNA consensus string into a MEME-format PWM motif.

Usage:
  python iupac2meme.py WTTTAYRTTTW > ACS_consensus.meme
  python iupac2meme.py WWWWTTTAYRTTTWGTT --name EACS_17bp_consensus > EACS_consensus.meme
  python iupac2meme.py ANWWAAAT --name B2_consensus > B2_consensus.meme
"""

import argparse

IUPAC_TO_PROBS = {
    "A": (1, 0, 0, 0),
    "C": (0, 1, 0, 0),
    "G": (0, 0, 1, 0),
    "T": (0, 0, 0, 1),
    "R": (0.5, 0, 0.5, 0),          # A/G
    "Y": (0, 0.5, 0, 0.5),          # C/T
    "S": (0, 0.5, 0.5, 0),          # C/G
    "W": (0.5, 0, 0, 0.5),          # A/T
    "K": (0, 0, 0.5, 0.5),          # G/T
    "M": (0.5, 0.5, 0, 0),          # A/C
    "B": (0, 1/3, 1/3, 1/3),        # C/G/T
    "D": (1/3, 0, 1/3, 1/3),        # A/G/T
    "H": (1/3, 1/3, 0, 1/3),        # A/C/T
    "V": (1/3, 1/3, 1/3, 0),        # A/C/G
    "N": (0.25, 0.25, 0.25, 0.25),  # A/C/G/T
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("iupac", help="IUPAC DNA consensus string, e.g. WTTTAYRTTTW")
    ap.add_argument("--name", default=None, help="Motif name (default: consensus string)")
    ap.add_argument("--bg", default="0.25 0.25 0.25 0.25",
                    help='Background frequencies as "A C G T" (default: uniform)')
    ap.add_argument("--strands", default="+ -", help='Strands line (default: "+ -")')
    args = ap.parse_args()

    s = args.iupac.strip().upper().replace("U", "T")
    name = args.name or s

    for ch in s:
        if ch not in IUPAC_TO_PROBS:
            raise SystemExit(f"Unsupported character '{ch}'. Allowed: {''.join(sorted(IUPAC_TO_PROBS))}")

    # Parse background
    bg_vals = [float(x) for x in args.bg.split()]
    if len(bg_vals) != 4 or any(v < 0 for v in bg_vals) or abs(sum(bg_vals) - 1.0) > 1e-6:
        raise SystemExit('Background must be 4 non-negative numbers summing to 1, e.g. "0.25 0.25 0.25 0.25"')
    bgA, bgC, bgG, bgT = bg_vals

    print("MEME version 4\n")
    print("ALPHABET= ACGT\n")
    print(f"strands: {args.strands}\n")
    print("Background letter frequencies")
    print(f"A {bgA:.6g} C {bgC:.6g} G {bgG:.6g} T {bgT:.6g}\n")
    print(f"MOTIF {name}")
    print(f"letter-probability matrix: alength= 4 w= {len(s)}")
    for ch in s:
        a, c, g, t = IUPAC_TO_PROBS[ch]
        #print(f"{a:.6g}\t{c:.6g}\t{g:.6g}\t{t:.6g}\t# {ch}")
        print(f"{a:.6g}\t{c:.6g}\t{g:.6g}\t{t:.6g}")

if __name__ == "__main__":
    main()
