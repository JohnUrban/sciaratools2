#!/usr/bin/env python3
"""
scan_protease_inhibitor_cys_motifs.py

Search protein FASTA sequences for cysteine-rich protease-inhibitor motifs
based on Fig. 1 from Dong, Xia, Zhao 2023:
"Antimicrobial components in the cocoon silk of silkworm, Bombyx mori"

Motifs:
  Kazal, Kunitz, TIL, WAP, pacifastin, I_68, I_83

Important:
  P1 is the protease-reactive-site position, NOT literal proline.
  This script captures whatever amino acid occurs at P1.

Outputs:
  1. per-protein summary TSV
  2. per-hit detail TSV
  3. optional hit FASTA

Coordinates are 1-based inclusive.

Example Usage:
python scan_protease_inhibitor_cys_motifs.py \
  -i sciara_proteins.faa \
  -o cys_inhibitor_scan/sciara_cys_inhibitor_scan \
  --hits-fasta \
  --include-negative-hit-rows

Produces:
cys_inhibitor_scan/sciara_cys_inhibitor_scan.protein_summary.tsv
cys_inhibitor_scan/sciara_cys_inhibitor_scan.motif_hits.tsv
cys_inhibitor_scan/sciara_cys_inhibitor_scan.motif_hits.faa

The summary table gives the per-protein classification:
protein_id  protein_length  has_motif  n_hits  motifs_detected  hit_names_detected  p1_residues_by_motif
geneA       183             YES        1       Kazal            Kazal-R             Kazal:R
geneB       220             NO         0       NA               NA                  NA

The hit table gives the detailed motif placement:
protein_id  has_motif  motif  hit_name  p1_residue  start_1based  end_1based  landmark_positions_1based
geneA       YES        Kazal  Kazal-R   R           24            91          C1=24;C2=30;P1=32;C3=40;...


Motifs grabbed from:
Antimicrobial components in the cocoon silk of silkworm, Bombyx mori. 2023.

"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional


@dataclass(frozen=True)
class Motif:
    name: str
    tokens: List[Tuple[str, str]]
    description: str = ""


# Token syntax:
#   ("C1", "C")         = conserved cysteine, literal C
#   ("c1", "C")         = relatively conserved cysteine, still literal C
#   ("P1", "[A-Z]")     = capture any amino acid at P1
#   ("x1", "[A-Z]{m,n}") = spacer
#
# The labels are used for output coordinates.
# Lowercase c labels preserve the figure's "relatively conserved cysteine" distinction.
MOTIFS: List[Motif] = [
    Motif(
        name="Kazal",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{1,4}"),
            ("C2", "C"),
            ("x2", "[A-Z]"),
            ("P1", "[A-Z]"),
            ("x3", "[A-Z]{5,10}"),
            ("C3", "C"),
            ("x4", "[A-Z]{10,14}"),
            ("C4", "C"),
            ("x5", "[A-Z]{6}"),
            ("C5", "C"),
            ("x6", "[A-Z]{10,20}"),
            ("C6", "C"),
        ],
        description="Kazal-like cysteine/P1 scaffold",
    ),
    Motif(
        name="Kunitz",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{8,9}"),
            ("c1", "C"),
            ("P1", "[A-Z]"),
            ("x2", "[A-Z]{14,18}"),
            ("C2", "C"),
            ("x3", "[A-Z]{7,11}"),
            ("c2", "C"),
            ("x4", "[A-Z]{11,12}"),
            ("C3", "C"),
            ("x5", "[A-Z]{3}"),
            ("C4", "C"),
        ],
        description="Kunitz-like cysteine/P1 scaffold",
    ),
    Motif(
        name="TIL",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{8,12}"),
            ("C2", "C"),
            ("x2", "[A-Z]{2,4}"),
            ("C3", "C"),
            ("x3", "[A-Z]{3,4}"),
            ("C4", "C"),
            ("x4", "[A-Z]{7,10}"),
            ("C5", "C"),
            ("x5", "[A-Z]"),
            ("P1", "[A-Z]"),
            ("x6", "[A-Z]{1,8}"),
            ("c1", "C"),
            ("x7", "[A-Z]{3,5}"),
            ("C6", "C"),
            ("x8", "[A-Z]"),
            ("C7", "C"),
            ("x9", "[A-Z]{11,13}"),
            ("C8", "C"),
            ("x10", "[A-Z]{5}"),
            ("C9", "C"),
        ],
        description="TIL-like cysteine/P1 scaffold",
    ),
    Motif(
        name="WAP",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{9,11}"),
            ("C2", "C"),
            ("P1", "[A-Z]"),
            ("x2", "[A-Z]{2}"),
            ("C3", "C"),
            ("x3", "[A-Z]{5}"),
            ("C4", "C"),
            ("x4", "[A-Z]{5}"),
            ("C5", "C"),
            ("C6", "C"),
            ("x5", "[A-Z]{3,4}"),
            ("C7", "C"),
            ("x6", "[A-Z]{3,4}"),
            ("C8", "C"),
        ],
        description="WAP-like cysteine/P1 scaffold",
    ),
    Motif(
        name="pacifastin",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{9,12}"),
            ("C2", "C"),
            ("x2", "[A-Z]{2}"),
            ("C3", "C"),
            ("x3", "[A-Z]"),
            ("C4", "C"),
            ("x4", "[A-Z]{5}"),
            ("C5", "C"),
            ("C6", "C"),
            ("x5", "[A-Z]{7,10}"),
            ("C7", "C"),
            ("x6", "[A-Z]"),
            ("P1", "[A-Z]"),
            ("x7", "[A-Z]{2,4}"),
            ("C8", "C"),
        ],
        description="Pacifastin-like cysteine/P1 scaffold",
    ),
    Motif(
        name="I_68",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{6}"),
            ("C2", "C"),
            ("x2", "[A-Z]{5}"),
            ("C3", "C"),
            ("x3", "[A-Z]{10,11}"),
            ("C4", "C"),
            ("x4", "[A-Z]{3,9}"),
            ("C5", "C"),
            ("C6", "C"),
            ("x5", "[A-Z]{6,9}"),
            ("C7", "C"),
            ("x6", "[A-Z]{6}"),
            ("C8", "C"),
            ("x7", "[A-Z]{3,6}"),
            ("C9", "C"),
            ("x8", "[A-Z]{9,11}"),
            ("C10", "C"),
            ("x9", "[A-Z]{5}"),
            ("C11", "C"),
            ("C12", "C"),
            ("x10", "[A-Z]{2}"),
            ("P1", "[A-Z]"),
        ],
        description="I_68-like cysteine/P1 scaffold",
    ),
    Motif(
        name="I_83",
        tokens=[
            ("C1", "C"),
            ("x1", "[A-Z]{4}"),
            ("C2", "C"),
            ("x2", "[A-Z]{4}"),
            ("C3", "C"),
            ("x3", "[A-Z]{6,8}"),
            ("C4", "C"),
            ("x4", "[A-Z]{13}"),
            ("C5", "C"),
            ("x5", "[A-Z]"),
            ("C6", "C"),
            ("C7", "C"),
            ("x6", "[A-Z]{2}"),
            ("C8", "C"),
            ("x7", "[A-Z]{9}"),
            ("C9", "C"),
            ("P1", "[A-Z]"),
            ("x8", "[A-Z]{13}"),
            ("C10", "C"),
            ("x9", "[A-Z]{5}"),
            ("C11", "C"),
            ("x10", "[A-Z]{4}"),
            ("C12", "C"),
        ],
        description="I_83-like cysteine/P1 scaffold",
    ),
]


def parse_fasta(path: Path) -> Iterable[Tuple[str, str, str]]:
    """Yield (record_id, description, sequence)."""
    header: Optional[str] = None
    seq_parts: List[str] = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    desc = header[1:].strip()
                    rec_id = desc.split()[0]
                    yield rec_id, desc, "".join(seq_parts).upper()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if header is not None:
        desc = header[1:].strip()
        rec_id = desc.split()[0]
        yield rec_id, desc, "".join(seq_parts).upper()


def build_regex(motif: Motif) -> re.Pattern:
    """Compile a motif regex with named capture groups for all tokens."""
    pieces = []
    for idx, (label, pattern) in enumerate(motif.tokens):
        group_name = f"g{idx}_{label.replace('-', '_')}"
        pieces.append(f"(?P<{group_name}>{pattern})")
    return re.compile("".join(pieces))


def token_group_names(motif: Motif) -> List[Tuple[str, str]]:
    """Return (label, regex group name) for motif tokens."""
    out = []
    for idx, (label, _pattern) in enumerate(motif.tokens):
        group_name = f"g{idx}_{label.replace('-', '_')}"
        out.append((label, group_name))
    return out


def find_overlapping(pattern: re.Pattern, seq: str) -> Iterable[re.Match]:
    """
    Find overlapping regex matches by trying each start position.

    This is slower than re.finditer, but protein FASTAs are usually small enough
    and this avoids missing overlapping motif occurrences.
    """
    for start in range(len(seq)):
        m = pattern.match(seq, start)
        if m is not None:
            yield m


def landmark_positions(m: re.Match, motif: Motif) -> Dict[str, int]:
    """
    Return 1-based positions for C/c/P1 landmarks only.
    If a motif had repeated labels this would need changing, but labels are unique here.
    """
    positions = {}
    for label, group_name in token_group_names(motif):
        if label.startswith("C") or label.startswith("c") or label == "P1":
            positions[label] = m.start(group_name) + 1
    return positions


def spacer_lengths(m: re.Match, motif: Motif) -> Dict[str, int]:
    """Return observed lengths for spacer groups."""
    lengths = {}
    for label, group_name in token_group_names(motif):
        if label.startswith("x"):
            lengths[label] = len(m.group(group_name))
    return lengths


def p1_residue(m: re.Match, motif: Motif) -> str:
    for label, group_name in token_group_names(motif):
        if label == "P1":
            return m.group(group_name)
    return "NA"


def format_positions(pos: Dict[str, int]) -> str:
    return ";".join(f"{k}={v}" for k, v in pos.items())


def format_lengths(lengths: Dict[str, int]) -> str:
    return ";".join(f"{k}={v}" for k, v in lengths.items())


def hit_fasta_header(
    protein_id: str,
    motif_name: str,
    p1: str,
    start: int,
    end: int,
    positions: Dict[str, int],
) -> str:
    pos_part = ":".join(f"{k}{v}" for k, v in positions.items())
    return f"{protein_id}|{motif_name}-{p1}|{start}-{end}|{pos_part}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan protein FASTA for cysteine-rich protease inhibitor motifs."
    )
    parser.add_argument("-i", "--input", required=True, type=Path, help="Protein FASTA")
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        type=Path,
        help="Output prefix, e.g. results/cys_motif_scan",
    )
    parser.add_argument(
        "--include-negative-hit-rows",
        action="store_true",
        help="Include one negative row in hit table for proteins with no motif hits.",
    )
    parser.add_argument(
        "--hits-fasta",
        action="store_true",
        help="Write FASTA of matched motif peptides.",
    )
    args = parser.parse_args()

    out_prefix: Path = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_path = out_prefix.with_suffix(".protein_summary.tsv")
    hits_path = out_prefix.with_suffix(".motif_hits.tsv")
    fasta_path = out_prefix.with_suffix(".motif_hits.faa")

    compiled = [(motif, build_regex(motif)) for motif in MOTIFS]

    summary_rows: List[List[str]] = []
    hit_rows: List[List[str]] = []
    fasta_records: List[Tuple[str, str]] = []

    for protein_id, desc, seq in parse_fasta(args.input):
        # Remove common non-AA symbols but keep X if present.
        clean_seq = re.sub(r"[^A-Z]", "", seq.upper())

        protein_hits = []

        for motif, pattern in compiled:
            for m in find_overlapping(pattern, clean_seq):
                start = m.start() + 1
                end = m.end()
                match_seq = m.group(0)
                p1 = p1_residue(m, motif)
                positions = landmark_positions(m, motif)
                lengths = spacer_lengths(m, motif)
                hit_name = f"{motif.name}-{p1}"

                protein_hits.append(
                    {
                        "protein_id": protein_id,
                        "protein_description": desc,
                        "motif": motif.name,
                        "hit_name": hit_name,
                        "p1_residue": p1,
                        "start": start,
                        "end": end,
                        "length": len(match_seq),
                        "landmark_positions": format_positions(positions),
                        "spacer_lengths": format_lengths(lengths),
                        "matched_sequence": match_seq,
                    }
                )

                if args.hits_fasta:
                    header = hit_fasta_header(
                        protein_id=protein_id,
                        motif_name=motif.name,
                        p1=p1,
                        start=start,
                        end=end,
                        positions=positions,
                    )
                    fasta_records.append((header, match_seq))

        if protein_hits:
            motif_names = sorted({h["motif"] for h in protein_hits})
            hit_names = sorted({h["hit_name"] for h in protein_hits})
            p1s_by_motif = {}
            for h in protein_hits:
                p1s_by_motif.setdefault(h["motif"], set()).add(h["p1_residue"])
            p1_summary = ";".join(
                f"{motif}:{','.join(sorted(p1s))}"
                for motif, p1s in sorted(p1s_by_motif.items())
            )

            summary_rows.append(
                [
                    protein_id,
                    desc,
                    str(len(clean_seq)),
                    "YES",
                    str(len(protein_hits)),
                    ",".join(motif_names),
                    ",".join(hit_names),
                    p1_summary,
                ]
            )

            for h in protein_hits:
                hit_rows.append(
                    [
                        h["protein_id"],
                        h["protein_description"],
                        "YES",
                        h["motif"],
                        h["hit_name"],
                        h["p1_residue"],
                        str(h["start"]),
                        str(h["end"]),
                        str(h["length"]),
                        h["landmark_positions"],
                        h["spacer_lengths"],
                        h["matched_sequence"],
                    ]
                )
        else:
            summary_rows.append(
                [
                    protein_id,
                    desc,
                    str(len(clean_seq)),
                    "NO",
                    "0",
                    "NA",
                    "NA",
                    "NA",
                ]
            )
            if args.include_negative_hit_rows:
                hit_rows.append(
                    [
                        protein_id,
                        desc,
                        "NO",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                    ]
                )

    with summary_path.open("w") as out:
        out.write(
            "\t".join(
                [
                    "protein_id",
                    "protein_description",
                    "protein_length",
                    "has_motif",
                    "n_hits",
                    "motifs_detected",
                    "hit_names_detected",
                    "p1_residues_by_motif",
                ]
            )
            + "\n"
        )
        for row in summary_rows:
            out.write("\t".join(row) + "\n")

    with hits_path.open("w") as out:
        out.write(
            "\t".join(
                [
                    "protein_id",
                    "protein_description",
                    "has_motif",
                    "motif",
                    "hit_name",
                    "p1_residue",
                    "start_1based",
                    "end_1based",
                    "hit_length",
                    "landmark_positions_1based",
                    "spacer_lengths",
                    "matched_sequence",
                ]
            )
            + "\n"
        )
        for row in hit_rows:
            out.write("\t".join(row) + "\n")

    if args.hits_fasta:
        with fasta_path.open("w") as out:
            for header, seq in fasta_records:
                out.write(f">{header}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + "\n")

    print(f"Wrote: {summary_path}")
    print(f"Wrote: {hits_path}")
    if args.hits_fasta:
        print(f"Wrote: {fasta_path}")


if __name__ == "__main__":
    main()
