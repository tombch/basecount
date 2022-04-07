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
20             1          0      1      0      0      0             0.0     100.0   0.0     0.0     0.0          0.0      0.0
21             1          1      0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
22             1          0      0      1      0      0             0.0     0.0     100.0   0.0     0.0          0.0      0.0
23             3          0      1      2      0      0             0.0     33.333  66.667  0.0     0.0          0.395    0.132
24             1334       0      0      0      1334   0             0.0     0.0     0.0     100.0   0.0          0.0      0.0
25             1396       1396   0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
26             1411       1411   0      0      0      0             100.0   0.0     0.0     0.0     0.0          0.0      0.0
27             1429       0      1427   0      2      0             0.0     99.86   0.0     0.14    0.0          0.007    0.0
28             1465       1453   0      1      0      11            99.181  0.0     0.068   0.0     0.751        0.031    0.0
29             1472       1471   0      1      0      0             99.932  0.0     0.068   0.0     0.0          0.004    0.0
30             1478       1464   1      13     0      0             99.053  0.068   0.88    0.0     0.0          0.035    0.0
:
```

#### Summary stats:
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME --summarise
```

Example output:
```
ref_name    first_ref_pos  last_ref_pos  avg_num_reads  avg_num_skip_del  avg_pc_skip_del  avg_entropy  avg_entropy_per_read
XXXXXXXXXX  20             29855         1623.174       34.149            2.04             0.108        0.001
```
