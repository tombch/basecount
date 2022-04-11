# `basecount`
## Usage
Create and activate the conda environment:
```
$ conda env create -f environment.yml
$ conda activate basecount
```

#### Per-position stats:
**All references**
```
$ python basecount.py BAM_FILE
```

**Selected references**
```
$ python basecount.py BAM_FILE --references REF_NAME_1 REF_NAME_2 ... 
```

Example output:
```
reference_name  reference_position  num_reads  num_a  num_c  num_g  num_t  num_deletions  pc_a    pc_c    pc_g    pc_t    pc_deletions  entropy  entropy_per_read
XXXXXXXXXX      0                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      1                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      2                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      3                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      4                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      5                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      6                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      7                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      8                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      9                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      10                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      11                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      12                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      13                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      14                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      15                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      16                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      17                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      18                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      19                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
XXXXXXXXXX      20                  1          0      1      0      0      0              0.0     100.0   0.0     0.0     0.0           0.0      0.0
XXXXXXXXXX      21                  1          1      0      0      0      0              100.0   0.0     0.0     0.0     0.0           0.0      0.0
XXXXXXXXXX      22                  1          0      0      1      0      0              0.0     0.0     100.0   0.0     0.0           0.0      0.0
XXXXXXXXXX      23                  3          0      1      2      0      0              0.0     33.333  66.667  0.0     0.0           0.395    0.132
XXXXXXXXXX      24                  1334       0      0      0      1334   0              0.0     0.0     0.0     100.0   0.0           0.0      0.0
XXXXXXXXXX      25                  1396       1396   0      0      0      0              100.0   0.0     0.0     0.0     0.0           0.0      0.0
XXXXXXXXXX      26                  1411       1411   0      0      0      0              100.0   0.0     0.0     0.0     0.0           0.0      0.0
XXXXXXXXXX      27                  1429       0      1427   0      2      0              0.0     99.86   0.0     0.14    0.0           0.007    0.0
XXXXXXXXXX      28                  1465       1453   0      1      0      11             99.181  0.0     0.068   0.0     0.751         0.031    0.0
XXXXXXXXXX      29                  1472       1471   0      1      0      0              99.932  0.0     0.068   0.0     0.0           0.004    0.0
:
```

#### Summary stats:
```
$ python basecount.py BAM_FILE --summarise
```

Example output:
```
ref_name                             XXXXXXXXXX
ref_length                           29903
num_reads                            123901
avg_coverage                         1649.805
pc_ref_coverage                      99.057
num_pos_no_coverage                  282
avg_num_deletions                    34.822
avg_pc_deletions                     2.011
avg_entropy                          0.116
avg_entropy_per_read                 0.01
mean_entropy_tile_vector             -
median_entropy_tile_vector           -
mean_entropy_per_read_tile_vector    -
median_entropy_per_read_tile_vector  -
```

#### Summary stats (with BED file):
```
$ python basecount.py BAM_PATH --bed BED_PATH --summarise
```

Example output:
```
ref_name                             XXXXXXXXXX
ref_length                           29903
num_reads                            123901
avg_coverage                         1649.805
pc_ref_coverage                      99.057
num_pos_no_coverage                  282
avg_num_deletions                    34.822
avg_pc_deletions                     2.011
avg_entropy                          0.116
avg_entropy_per_read                 0.01
mean_entropy_tile_vector             0.107, 0.116, 0.127, 0.118, 0.113, 0.111, 0.11, 0.077, 0.12, 0.071, 0.105, 0.119, 0.088, 0.147, 0.116, 0.109,..
median_entropy_tile_vector           0.094, 0.104, 0.112, 0.102, 0.108, 0.101, 0.104, 0.0, 0.113, 0.013, 0.093, 0.111, 0.022, 0.143, 0.107, 0.098,..
mean_entropy_per_read_tile_vector    0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006, 0.0, 0.003, 0.0, 0.0, 0.006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0..
median_entropy_per_read_tile_vector  0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0...
```
