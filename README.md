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

COPYING
---------
Copyright (c) 2021 jnjcc, [Yste.org](http://www.yste.org)

This program is licensed under the terms of the MIT license. See [COPYING](COPYING)
for more details.
