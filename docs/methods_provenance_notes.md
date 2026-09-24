# Methods provenance notes

## HiFi extractor implementation

The authoritative development-archive G6 closeout committed the frozen high-accuracy truth rules and the executed HiFi generator together. One implementation-level sentence in frozen rule H05 is stale: it describes exact 11-mer seeding plus banded semi-global anchor matching.

The tracked generator `analysis/frozen_source/g6/g6_hifi_extract.py`, which generated the frozen HiFi evidence used by the final result tables, implements:
- 25-bp left/right catalogue anchors;
- exact 16-mer seed indexing;
- direct Hamming verification with at most two mismatches per anchor;
- both orientations;
- tract spans from 0 to 150 bp;
- cross-locus arbitration before evidence emission.

The contemporaneous G6 high-accuracy rescue report in the development archive explicitly describes the executed extractor as an own implementation using exact k-mer seeding with Hamming anchors.

Publication policy: report the tracked generator as the executed implementation; preserve the H05 prose unchanged as historical provenance; do not re-extract reads or change frozen truth thresholds during manuscript preparation. H06-H13 truth-status and comparison criteria are unaffected.

## G7 follow-up bounded execution

The frozen PRJNA678459 follow-up is intentionally bounded and disease-blind. Caller extraction was capped at 500,000 read fragments per sample. Screening-structure sketches used canonical 21-mers, bottom-1,000 deterministic crc32 hashes, and the first 30,000 read pairs per sample. These caps are implementation scope limits and must be disclosed wherever the corresponding pilot or structure outputs are described.
