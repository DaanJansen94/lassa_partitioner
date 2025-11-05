# LassaPartitioner

A bioinformatics tool specifically designed for creating APOBEC3 and non-APOBEC3 partitions from Lassa virus sequence alignments. This tool detects APOBEC3 editing signatures (GA→AA and TC→TT conversions) that are common in Lassa virus genomes.

## Background

APOBEC3 deaminases can edit viral genomes, producing characteristic mutation patterns. In Lassa virus alignments, we can observe GA→AA and TC→TT conversions that indicate APOBEC3 activity. When performing evolutionary analysis, it's important to separate these APOBEC3-edited sites from other mutations to avoid bias from sequencing errors or bioinformatics artifacts.

This tool identifies positions with evidence of APOBEC3 conversion activity (variation between unconverted and converted states) and creates two complementary partitions for phylogenetic analysis.

## Key Features

- **Detects APOBEC3 conversion activity**: Only includes sites with evidence of GA→AA or TC→TT conversions (variation between sequences)
- **Lassa-specific optimization**: Optimized for detecting APOBEC3 signatures in Lassa virus alignments
- **Statistics output**: Automatically generates statistics files showing the percentage of APOBEC3 and non-APOBEC3 sites
- **Dual partitions**: Creates both APOBEC3 and non-APOBEC3 partitions for comprehensive analysis

## Installation

### Option 1: From Source Code

1. Ensure you have Python 3.6 or higher installed:
   ```bash
   python3 --version
   ```

2. Install Biopython:
   ```bash
   pip install biopython
   ```

3. Clone and install LassaPartitioner:
   ```bash
   git clone https://github.com/DaanJansen94/lassa_partitioner
   cd lassa_partitioner
   pip install -e .
   ```

### Option 2: Using pip

```bash
pip install lassa-partitioner
```

## Usage

### Basic Command Structure

```bash
LassaPartitioner -partition <type> -i <input_file> -o <output_file>
```

### Options

- `-partition`: Type of partition to create
  - `apobec`: Create only APOBEC3 partition
  - `non-apobec`: Create only non-APOBEC3 partition
  - `both`: Create both partitions (recommended)

- `-i, --input`: Input FASTA alignment file (Lassa virus sequences)

- `-o, --output`: Output file prefix

### Examples

#### 1. Create both partitions (recommended):
```bash
LassaPartitioner -partition both -i lassa_alignment.fasta -o lassa_partitioned
```

This creates:
- `lassa_partitioned_apobec.fasta`: APOBEC3 partition
- `lassa_partitioned_non-apobec.fasta`: Non-APOBEC3 partition
- `lassa_partitioned.stats.txt`: Statistics file

#### 2. Create only APOBEC3 partition:
```bash
LassaPartitioner -partition apobec -i lassa_alignment.fasta -o lassa_apobec
```

#### 3. Create only non-APOBEC3 partition:
```bash
LassaPartitioner -partition non-apobec -i lassa_alignment.fasta -o lassa_nonapobec
```

## Output Files

### When using `-partition both`:

1. **`{output}_apobec.fasta`**: APOBEC3 partition
   - Contains only positions with evidence of APOBEC3 conversion activity (GA→AA or TC→TT variation)
   - All other positions are masked as 'N'
   - Suitable for phylogenetic analysis focusing on APOBEC3-edited sites

2. **`{output}_non-apobec.fasta`**: Non-APOBEC3 partition
   - Contains all positions except APOBEC3 sites
   - APOBEC3 target positions are masked as 'N'
   - Complementary to the APOBEC3 partition

3. **`{output}.stats.txt`**: Statistics file
   - Shows the number and percentage of non-N bases in each partition
   - Example output:
     ```
     Partition Statistics
     ====================
     
     APOBEC3 Partition:
       Non-N bases: 371,691/805,134
       Percentage: 46.17%
     
     Non-APOBEC3 Partition:
       Non-N bases: 355,600/805,134
       Percentage: 44.17%
     ```

## How It Works

### APOBEC3 Detection

The tool identifies APOBEC3 sites by detecting:

1. **GA→AA conversions**: Positions where some sequences have `GA` dinucleotides and others have `AA` (indicating G→A conversion)
2. **TC→TT conversions**: Positions where some sequences have `TC` dinucleotides and others have `TT` (indicating C→T conversion)

**Important**: The tool only includes positions with **evidence of conversion activity** (variation between sequences). Positions with only GA or only TC (without conversion evidence) are masked as N, as they don't show active APOBEC3 editing.

### Example: GA→AA Conversion

At a position where APOBEC3 editing occurred:
- Sequence 1: `GA` (unconverted)
- Sequence 2: `GA` (unconverted)
- Sequence 3: `AA` (converted - G was edited to A)
- Sequence 4: `GA` (unconverted)

This position would be included in the APOBEC3 partition because it shows variation indicating active editing.

## Use Cases

### Phylogenetic Analysis

Use the APOBEC3 partition for:
- Evolutionary analysis focusing on APOBEC3-edited sites
- Removing bias from sequencing errors
- Studying APOBEC3 editing patterns in Lassa virus

Use the non-APOBEC3 partition for:
- Complementary analysis excluding APOBEC3 sites
- Studying mutations from other sources (polymerase errors, etc.)

### Combining Partitions

You can use both partitions in separate phylogenetic analyses or combine them depending on your research question.

## Requirements

- Python 3.6 or higher
- Biopython >= 1.80, <= 1.90

## Citation

If you use LassaPartitioner in your research, please cite:

```
Jansen, D. (2025). LassaPartitioner: A bioinformatics tool for creating APOBEC3 and 
non-APOBEC3 partitions from Lassa virus sequence alignments. GitHub. 
https://github.com/DaanJansen94/lassa_partitioner
```

## License

This project is licensed under MIT License - see the [LICENSE](LICENSE) file for details.

## Support

If you encounter any problems or have questions, please open an issue on [GitHub](https://github.com/DaanJansen94/lassa_partitioner/issues).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

