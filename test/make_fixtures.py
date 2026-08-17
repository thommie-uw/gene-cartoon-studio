#!/usr/bin/env python3
"""Re-capture the genePred rows used by test_offline.py from the live UCSC API.

Run this if UCSC updates its annotation and you want the offline test to track
the new models:

    python test/make_fixtures.py > test/fixtures/latest.json

then paste the rows of interest into test_offline.py.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ucsc_gene_cartoon import UCSCClient  # noqa: E402

REGIONS = [
    ("TP53", "ncbiRefSeqSelect", "chr17", 7668000, 7688000),
    ("TP53-isoforms", "ncbiRefSeqCurated", "chr17", 7668000, 7688000),
    ("CFTR", "ncbiRefSeqSelect", "chr7", 117480000, 117670000),
    ("XIST", "ncbiRefSeqSelect", "chrX", 73815000, 73860000),
]

if __name__ == "__main__":
    c = UCSCClient("hg38", quiet=True)
    out = {}
    for label, track, chrom, start, end in REGIONS:
        txs = c.transcripts_in_region(track, chrom, start, end)
        out[label] = [t.__dict__ for t in txs]
        print(f"{label}: {len(txs)} transcript(s) from {track}", file=sys.stderr)
    json.dump(out, sys.stdout, indent=2, default=str)
