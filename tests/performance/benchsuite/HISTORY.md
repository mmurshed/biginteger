# BigMath vs GMP — historical comparison

Selected large-size results, one block per run, appended by `run_bigmath`.
Lower ratio is better (BigMath ms / GMP ms). Rows are the entries flagged
`historical` in the dataset manifest.

## run 1 — machine=Apple_M1_Max_arm64 profile=quick

| op | size | BigMath ms | GMP ms | BM/GMP | check |
|----|------|-----------:|-------:|-------:|:-----:|
| mul | 10000x10000 digits | 0.093 | 0.030 | 3.12x | FAIL |
| div | 40000x10000 digits | 1.186 | 0.232 | 5.10x | FAIL |

## run 2 — machine=Apple_M1_Max_arm64 profile=quick

| op | size | BigMath ms | GMP ms | BM/GMP | check |
|----|------|-----------:|-------:|-------:|:-----:|
| mul | 10000x10000 digits | 0.231 | 0.030 | 7.76x | ok |
| div | 40000x10000 digits | 1.545 | 0.232 | 6.65x | ok |

## run 1 — machine=Apple_M1_Max_arm64 profile=default

| op | size | BigMath ms | GMP ms | BM/GMP | check |
|----|------|-----------:|-------:|-------:|:-----:|
| mul | 1000000x1000000 digits | 10.682 | 8.716 | 1.23x | ok |
| div | 1000000x200000 digits | 34.660 | 10.367 | 3.34x | ok |

## run 2 — machine=Apple_M1_Max_arm64 profile=default — DISCARD (ambient load avg 6–8; all BigMath rows ~2.3× inflated incl. untouched mul)

| op | size | BigMath ms | GMP ms | BM/GMP | check |
|----|------|-----------:|-------:|-------:|:-----:|
| mul | 1000000x1000000 digits | 24.396 | 8.716 | 2.80x | ok |
| div | 1000000x200000 digits | 58.367 | 10.367 | 5.63x | ok |

## run 3 — machine=Apple_M1_Max_arm64 profile=default — DISCARD (stale cached run_bigmath binary from before PRs #82–#89; measured pre-session code)

| op | size | BigMath ms | GMP ms | BM/GMP | check |
|----|------|-----------:|-------:|-------:|:-----:|
| mul | 1000000x1000000 digits | 10.816 | 8.716 | 1.24x | ok |
| div | 1000000x200000 digits | 33.872 | 10.367 | 3.27x | ok |

## run 4 — machine=Apple_M1_Max_arm64 profile=default — post PRs #82–#89 (FastDiv shift-norm, D&C thresholds, wraparound Newton family, balanced band 4/3 + quotient-sized division)

| op | size | BigMath ms | GMP ms | BM/GMP | check |
|----|------|-----------:|-------:|-------:|:-----:|
| mul | 1000000x1000000 digits | 10.641 | 8.716 | 1.22x | ok |
| div | 1000000x200000 digits | 22.706 | 10.367 | 2.19x | ok |
