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
reference_position  num_reads  num_a  num_c  num_g  num_t  num_deletions  pc_a    pc_c    pc_g    pc_t    pc_deletions  entropy  entropy_per_read
0                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
1                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
2                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
3                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
4                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
5                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
6                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
7                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
8                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
9                   0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
10                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
11                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
12                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
13                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
14                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
15                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
16                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
17                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
18                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
19                  0          0      0      0      0      0              -1      -1      -1      -1      -1            1        1
20                  1          0      1      0      0      0              0.0     100.0   0.0     0.0     0.0           0.0      0.0
21                  1          1      0      0      0      0              100.0   0.0     0.0     0.0     0.0           0.0      0.0
22                  1          0      0      1      0      0              0.0     0.0     100.0   0.0     0.0           0.0      0.0
23                  3          0      1      2      0      0              0.0     33.333  66.667  0.0     0.0           0.395    0.132
24                  1334       0      0      0      1334   0              0.0     0.0     0.0     100.0   0.0           0.0      0.0
25                  1396       1396   0      0      0      0              100.0   0.0     0.0     0.0     0.0           0.0      0.0
26                  1411       1411   0      0      0      0              100.0   0.0     0.0     0.0     0.0           0.0      0.0
27                  1429       0      1427   0      2      0              0.0     99.86   0.0     0.14    0.0           0.007    0.0
28                  1465       1453   0      1      0      11             99.181  0.0     0.068   0.0     0.751         0.031    0.0
29                  1472       1471   0      1      0      0              99.932  0.0     0.068   0.0     0.0           0.004    0.0
:
```

#### Summary stats:
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME --summarise
```

Example output:
```
ref_name                             XXXXXXXXXX
ref_length                           29903
num_no_coverage                      282
pc_no_coverage                       0.943
avg_num_reads                        1649.805
avg_num_deletions                    34.822
avg_pc_deletions                     2.011
avg_entropy                          0.116
avg_entropy_per_read                 0.01
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
ref_length                           29903
num_no_coverage                      282
pc_no_coverage                       0.943
avg_num_reads                        1649.805
avg_num_deletions                    34.822
avg_pc_deletions                     2.011
avg_entropy                          0.116
avg_entropy_per_read                 0.01
median_entropy_tile_vector           0.094,0.104,0.112,0.102,0.108,0.101,0.104,0.0,0.113,0.013,0.093,0.111,0.022,0.143,0.107,0.098,..
median_entropy_per_read_tile_vector  0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,..
```
