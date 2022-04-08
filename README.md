# `basecount`
## Usage
Create and activate the conda environment:
```
$ conda env create -f environment.yml
$ conda activate basecount
```

#### Per-position stats:
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME
```

Example output:
```
reference_pos  num_reads  num_a  num_c  num_g  num_t  num_skip_del  pc_a    pc_c    pc_g    pc_t    pc_skip_del  entropy  entropy_per_read
24             387        0      0      0      387    0             0.0     0.0     0.0     100.0   0.0          0.0      0.0
25             407        407    0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
26             409        409    0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
27             413        0      412    0      1      0             0.0     99.758  0.0     0.242   0.0          0.011    0.0
28             419        416    0      0      0      3             99.284  0.0     0.0     0.0     0.716        0.026    0.0
29             421        421    0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
30             423        421    0      2      0      0             99.527  0.0     0.473   0.0     0.0          0.019    0.0
31             425        0      423    1      0      1             0.0     99.529  0.235   0.0     0.235        0.021    0.0
32             426        0      425    0      1      0             0.0     99.765  0.0     0.235   0.0          0.01     0.0
33             426        426    0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
34             426        426    0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
:
```

#### Summary stats:
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME --summarise
```

Example output:
```
ref_name                             XXXXXXXXXX
first_ref_pos                        24
last_ref_pos                         29853
avg_num_reads                        210.74
avg_num_skip_del                     4.52
avg_pc_skip_del                      1.89
avg_entropy                          0.085
avg_entropy_per_read                 0.003
median_entropy_tile_vector           -
median_entropy_per_read_tile_vector  -
```

#### Summary stats (with BED file):
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME --bed BED_PATH --summarise
```

Example output:
```
ref_name                             XXXXXXXXXX
first_ref_pos                        24
last_ref_pos                         29853
avg_num_reads                        210.74
avg_num_skip_del                     4.52
avg_pc_skip_del                      1.89
avg_entropy                          0.085
avg_entropy_per_read                 0.003
median_entropy_tile_vector           0.0,0.082,0.113,0.102,0.0,0.095,0.071,0.081,0.0,0.101,0.071,0.099,0.067,0.088,0.09,0.098,0.104,0.0,0.0,..
median_entropy_per_read_tile_vector  0.0,0.001,0.003,0.0,0.0,0.0,0.002,0.002,0.0,0.001,0.001,0.001,0.001,0.001,0.001,0.0,0.004,0.0,0.0,0.002..
```
