# Quick Start Guide for LassaPartitioner

## Installation

```bash
git clone https://github.com/DaanJansen94/lassa_partitioner
cd lassa_partitioner
pip install -e .
```

## Basic Usage

### Create both partitions (recommended):

```bash
LassaPartitioner -partition both -i your_lassa_alignment.fasta -o output_prefix
```

This creates:
- `output_prefix_apobec.fasta` - APOBEC3 partition
- `output_prefix_non-apobec.fasta` - Non-APOBEC3 partition  
- `output_prefix.stats.txt` - Statistics file

### Example with Lassa L segment:

```bash
LassaPartitioner -partition both -i l_trimmed.fasta -o lassa_partitioned
```

## Output Explained

### APOBEC3 Partition
Contains only positions with evidence of APOBEC3 conversion activity:
- GA→AA conversions (where some sequences have GA, others have AA)
- TC→TT conversions (where some sequences have TC, others have TT)
- All other positions are masked as 'N'

### Non-APOBEC3 Partition
Contains all positions except APOBEC3 sites:
- APOBEC3 target positions are masked as 'N'
- Everything else is retained

### Statistics File
Shows the percentage of non-N bases in each partition, helping you understand how much data is retained.

## Key Feature

**Important**: LassaPartitioner only includes positions with **evidence of conversion activity**. If all sequences have GA (no AA) or all have TC (no TT), those positions are masked as N because they don't show active APOBEC3 editing.

This ensures you're only analyzing sites where APOBEC3 editing actually occurred, not just potential target sites.

## Need Help?

- Full documentation: See [README.md](README.md)
- Issues: https://github.com/DaanJansen94/lassa_partitioner/issues
