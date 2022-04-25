# `basecount`
## Usage
Create and activate the conda environment:
```
$ conda env create -f environment.yml
$ conda activate basecount
```

Compile the C++ code:
```
$ c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) count.cpp -o count$(python3-config --extension-suffix)
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
reference   position  coverage  num_a  num_c  num_g  num_t  num_n  num_ds  pc_a    pc_c    pc_g    pc_t    pc_n  pc_ds   entropy
MN908947.3  1         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  2         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  3         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  4         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  5         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  6         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  7         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  8         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  9         0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  10        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  11        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  12        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  13        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  14        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  15        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  16        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  17        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  18        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  19        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  20        0         0      0      0      0      0      0       -1      -1      -1      -1      -1    -1      1
MN908947.3  21        1         0      1      0      0      0      0       0.0     100.0   0.0     0.0     0.0   0.0     0.0
MN908947.3  22        1         1      0      0      0      0      0       100.0   0.0     0.0     0.0     0.0   0.0     0.0
MN908947.3  23        1         0      0      1      0      0      0       0.0     0.0     100.0   0.0     0.0   0.0     0.0
MN908947.3  24        3         0      1      2      0      0      0       0.0     33.333  66.667  0.0     0.0   0.0     0.395
MN908947.3  25        1334      0      0      0      1334   0      0       0.0     0.0     0.0     100.0   0.0   0.0     0.0
MN908947.3  26        1396      1396   0      0      0      0      0       100.0   0.0     0.0     0.0     0.0   0.0     0.0
MN908947.3  27        1411      1411   0      0      0      0      0       100.0   0.0     0.0     0.0     0.0   0.0     0.0
MN908947.3  28        1429      0      1427   0      2      0      0       0.0     99.86   0.0     0.14    0.0   0.0     0.007
MN908947.3  29        1465      1453   0      1      0      0      11      99.181  0.0     0.068   0.0     0.0   0.751   0.031
MN908947.3  30        1472      1471   0      1      0      0      0       99.932  0.0     0.068   0.0     0.0   0.0     0.004
MN908947.3  31        1478      1464   1      13     0      0      0       99.053  0.068   0.88    0.0     0.0   0.0     0.035
MN908947.3  32        1482      0      1479   0      2      0      1       0.0     99.798  0.0     0.135   0.0   0.067   0.01
MN908947.3  33        1483      1      1479   0      3      0      0       0.067   99.73   0.0     0.202   0.0   0.0     0.013
MN908947.3  34        1483      1482   0      1      0      0      0       99.933  0.0     0.067   0.0     0.0   0.0     0.003
MN908947.3  35        1483      1481   0      2      0      0      0       99.865  0.0     0.135   0.0     0.0   0.0     0.006
MN908947.3  36        1484      0      1480   0      2      0      2       0.0     99.73   0.0     0.135   0.0   0.135   0.013
MN908947.3  37        1485      0      1485   0      0      0      0       0.0     100.0   0.0     0.0     0.0   0.0     0.0
MN908947.3  38        1485      1485   0      0      0      0      0       100.0   0.0     0.0     0.0     0.0   0.0     0.0
MN908947.3  39        1485      1483   0      2      0      0      0       99.865  0.0     0.135   0.0     0.0   0.0     0.006
:
```
The table can also be displayed in long format (one base per row) by passing `--long` as an argument.

#### Summary stats:
```
$ python basecount.py BAM_FILE --summarise
```

Example output:
```
ref_name                    XXXXXXXXXX
ref_length                  29903
num_reads                   123901
avg_coverage                1649.805
pc_ref_coverage             99.057
num_pos_no_coverage         282
avg_num_deletions_skips     34.822
avg_pc_deletions_skips      2.011
avg_entropy                 0.116
mean_entropy_tile_vector    -
median_entropy_tile_vector  -
```

#### Summary stats (with BED file):
```
$ python basecount.py BAM_PATH --bed BED_PATH --summarise
```

Example output:
```
ref_name                    XXXXXXXXXX
ref_length                  29903
num_reads                   123901
avg_coverage                1649.805
pc_ref_coverage             99.057
num_pos_no_coverage         282
avg_num_deletions_skips     34.822
avg_pc_deletions_skips      2.011
avg_entropy                 0.116
mean_entropy_tile_vector    0.107, 0.116, 0.127, 0.118, 0.113, 0.111, 0.11, 0.077, 0.12, 0.071, 0.105, 0.119, 0.088, 0.147, 0.116, 0.109, 0.122, 0.12, ..
median_entropy_tile_vector  0.094, 0.104, 0.112, 0.102, 0.108, 0.101, 0.104, 0.0, 0.113, 0.013, 0.093, 0.111, 0.022, 0.143, 0.107, 0.098, 0.108, 0.106,..
```
