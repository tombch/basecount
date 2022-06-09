# `basecount`
### Setup
```
$ git clone https://github.com/tombch/basecount.git
$ cd basecount/
$ conda env create -f environment.yml
$ conda activate basecount
$ pip install .
```

### Usage

#### Per-position statistics:
```
$ basecount BAM_FILE
```

Example output:
```
reference   position  coverage  num_a  num_c  num_g  num_t  num_ds  pc_a    pc_c    pc_g    pc_t    pc_ds   entropy  secondary_entropy
XXXXXXXXXX  1         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  2         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  3         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  4         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  5         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  6         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  7         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  8         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  9         0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  10        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  11        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  12        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  13        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  14        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  15        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  16        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  17        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  18        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  19        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  20        0         0      0      0      0      0       -1      -1      -1      -1      -1      1        1
XXXXXXXXXX  21        1         0      1      0      0      0       0.0     100.0   0.0     0.0     0.0     0.0      1
XXXXXXXXXX  22        1         1      0      0      0      0       100.0   0.0     0.0     0.0     0.0     0.0      1
XXXXXXXXXX  23        1         0      0      1      0      0       0.0     0.0     100.0   0.0     0.0     0.0      1
XXXXXXXXXX  24        3         0      1      2      0      0       0.0     33.333  66.667  0.0     0.0     0.395    0.0
XXXXXXXXXX  25        1334      0      0      0      1334   0       0.0     0.0     0.0     100.0   0.0     0.0      1
XXXXXXXXXX  26        1396      1396   0      0      0      0       100.0   0.0     0.0     0.0     0.0     0.0      1
XXXXXXXXXX  27        1411      1411   0      0      0      0       100.0   0.0     0.0     0.0     0.0     0.0      1
XXXXXXXXXX  28        1429      0      1427   0      2      0       0.0     99.86   0.0     0.14    0.0     0.007    0.0
XXXXXXXXXX  29        1465      1453   0      1      0      11      99.181  0.0     0.068   0.0     0.751   0.031    0.207
XXXXXXXXXX  30        1472      1471   0      1      0      0       99.932  0.0     0.068   0.0     0.0     0.004    0.0
XXXXXXXXXX  31        1478      1464   1      13     0      0       99.053  0.068   0.88    0.0     0.0     0.035    0.186
XXXXXXXXXX  32        1482      0      1479   0      2      1       0.0     99.798  0.0     0.135   0.067   0.01     0.459
XXXXXXXXXX  33        1483      1      1479   0      3      0       0.067   99.73   0.0     0.202   0.0     0.013    0.406
XXXXXXXXXX  34        1483      1482   0      1      0      0       99.933  0.0     0.067   0.0     0.0     0.003    0.0
XXXXXXXXXX  35        1483      1481   0      2      0      0       99.865  0.0     0.135   0.0     0.0     0.006    0.0
XXXXXXXXXX  36        1484      0      1480   0      2      2       0.0     99.73   0.0     0.135   0.135   0.013    0.5
XXXXXXXXXX  37        1485      0      1485   0      0      0       0.0     100.0   0.0     0.0     0.0     0.0      1
XXXXXXXXXX  38        1485      1485   0      0      0      0       100.0   0.0     0.0     0.0     0.0     0.0      1
XXXXXXXXXX  39        1485      1483   0      2      0      0       99.865  0.0     0.135   0.0     0.0     0.006    0.0
:
```

#### Summary statistics:
```
$ basecount BAM_FILE --summarise
```

Example output:
```
reference_name         XXXXXXXXXX
reference_length       29903
num_reads              123901
pc_reference_coverage  99.057
avg_depth              1649.805
avg_entropy            0.116
```

#### Summary statistics (with BED file):
```
$ basecount BAM_FILE --summarise-with-bed BED_FILE
```

Example output:
```
reference_name                            XXXXXXXXXX
reference_length                          29903
num_reads                                 123901
pc_reference_coverage                     99.057
avg_depth                                 1649.805
avg_entropy                               0.116
mean_coverage_amplicon_vector             243.906, 1229.481, 1729.518, 2711.367, 1346.401, 3041.555, 1309.193, 143.864, 1494.996, 42..
median_coverage_amplicon_vector           138.0, 1219.0, 1725.0, 2582.0, 678.0, 2932.0, 1309.0, 13.0, 1496.0, 16.0, 1471.0, 1644.0, ..
mean_entropy_amplicon_vector              0.107, 0.116, 0.127, 0.118, 0.113, 0.111, 0.11, 0.077, 0.12, 0.071, 0.105, 0.119, 0.088, 0..
median_entropy_amplicon_vector            0.094, 0.104, 0.112, 0.102, 0.108, 0.101, 0.104, 0.0, 0.113, 0.013, 0.093, 0.111, 0.022, 0..
mean_secondary_entropy_amplicon_vector    0.437, 0.724, 0.702, 0.724, 0.69, 0.714, 0.709, 0.666, 0.713, 0.685, 0.735, 0.724, 0.666, ..
median_secondary_entropy_amplicon_vector  0.485, 0.739, 0.714, 0.749, 0.729, 0.746, 0.736, 1.0, 0.737, 0.899, 0.753, 0.736, 0.953, 0..
```

#### View all options:

```
$ basecount --help
```
