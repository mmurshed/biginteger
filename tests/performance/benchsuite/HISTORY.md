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
