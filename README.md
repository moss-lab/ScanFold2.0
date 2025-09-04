# ScanFold2.0: RNA Structure Scanning Tool

ScanFold2.0 is a tool for analyzing RNA structure through scanning window predictions, evaluating thermodynamic stability, and identifying statistically significant structural motifs. It uses machine learning models to calculate z-scores and performs competitive base pair resolution.

## Usage

```bash
python ScanFold2.0.py <input_fasta> [options]
```

## Parameters

### Basic Parameters

| Parameter | Description | Default Value | Notes |
|-----------|-------------|---------------|-------|
| `-w` | Window size in nucleotides | `120` | Length of the sliding window |
| `-s` | Step size in nucleotides | `1` | Distance between consecutive windows |
| `-t` | Folding temperature in °C | `37` | Temperature used for RNA folding |
| `-r` | Number of randomizations | `100` | Used for background shuffling |

### Structure Analysis Parameters

| Parameter | Description | Default Value | Notes |
|-----------|-------------|---------------|-------|
| `-f` | Z-score filter threshold | `-2` | Base pairs with Z-scores below this value are included |
| `-c` | Competition parameter | `1` | `1` = disallow competition (resolve competing pairs)<br>`0` = allow competition |
| `--global_refold` | Enable global refolding | `False` | Refolds full sequence using filtered base pairs |

### Output Parameters

| Parameter | Description | Default Value | Notes |
|-----------|-------------|---------------|-------|
| `--id` | Sequence identifier | `UserInput` | Name/ID used in output files |
| `--folder_name` | Output folder name | `None` | Uses header name or datetime if not specified |

### Algorithm Parameters

| Parameter | Description | Default Value | Notes |
|-----------|-------------|---------------|-------|
| `--algorithm` | RNA folding algorithm | `rnafold` | Options: `rnafold`, `rnastructure`, `mxfold` |

### WebServer Parameters

| Parameter | Description | Default Value | Notes |
|-----------|-------------|---------------|-------|
| `--webserver` | Path for webserver output | `None` | Only evaluates first sequence in FASTA if provided |
| `--logfile` | Path to write log file | `sys.stdout` | Output destination for logging information |
| `--loglevel` | Logging level | `INFO` | Standard Python logging levels |

## Output Files

ScanFold2.0 generates several output files:

- `.tsv` or `.csv` - Raw scanning data with z-scores
- `.ct` - Connectivity table files at various z-score thresholds
- `.dbn` - Dot-bracket notation files compatible with structure visualization tools
- `.bp` - Base pair information
- `.wig` - Wiggle files for genome browser visualization
- `.log` - Log files with analysis details

## Notes

- Lower Z-score values (more negative) indicate more stable and significant structures
- The competitive base pairing resolution (`-c 1`) is recommended for most analyses
- For visualizing structures, use the `.dbn` files with tools like Forna or the ScanFoldBothForInforna.py script

For detailed information about the algorithm, please see the original publications or contact the developers.
