# `basecount`
## Usage
Create and activate the conda environment:
```
$ conda env create -f environment.yml
$ conda activate basecount
```

Generate per-position stats:
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME
```

Generate summary stats:
```
$ python basecount.py --bam BAM_PATH --ref REFERENCE_NAME --summarise
```