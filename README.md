pymatch2
=====================
An extension of the [pymatch](https://github.com/benmiroglio/pymatch) package,
with faster Nearest Neighbor Matching algorithms.

Tests
---------
| Device           |   macOS | CPU                   | RAM | SSD   |
|------------------|---------|-----------------------|-----|-------|
| MacBook Pro 2017 | 10.13.6 | 2.3 GHz Intel Core i5 | 8GB | 256GB |

```python
## test performance of matching methods (in seconds)
python -m test.benchmark
```

| # test | # control | random |    min |    nnm | bsearch |
|--------|-----------|--------|--------|--------|---------|
|  10000 |    200000 |  15.85 | 287.85 |  27.19 |   20.50 |
|  15000 |    200000 |  23.27 | 421.99 |  29.65 |   33.65 |
|  20000 |    200000 |  32.99 | 581.42 |  32.70 |   41.42 |
|  50000 |    900000 | 285.76 |      - | 128.37 |  119.63 |
| 100000 |   1000000 | 659.55 |      - | 155.13 |  250.11 |

```python
python -m test.benchmark
```

| # test | # control |  random | min    |    nnm | bsearch |  merge |
|--------|-----------|---------|--------|--------|---------|--------|
|  10000 |    200000 |    5.63 | 193.79 |   9.70 |    7.82 |   9.91 |
|  15000 |    200000 |    8.62 | 293.11 |  10.80 |   11.70 |  10.73 |
|  20000 |    200000 |   11.45 | 393.03 |  11.67 |   15.68 |  11.63 |
|  50000 |    900000 |   90.26 |      - |  45.81 |   43.73 |  45.39 |
| 100000 |   1000000 |  198.50 |      - |  58.17 |   87.52 |  57.92 |
| 500000 |   5000000 | 4071.77 |      - | 293.15 |  491.03 | 290.33 |

COPYING
---------
Copyright (c) 2021 jnjcc, [Yste.org](http://www.yste.org)

This program is licensed under the terms of the MIT license. See [COPYING](COPYING)
for more details.
