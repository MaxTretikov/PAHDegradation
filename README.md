# Pseudomonas sp. JM2 PAH Degradation with bistable switch

Simulation of PAH degradation efficiency with an inducible vs bistable system.

## Installation and use

```
pip install -e git+https://github.com/MaxTretikov/PAHDegradation.git
```

Generate figure 1:
```
python -m pahdegradation
```

Generate figure 2:
```
python -m pahdegradation --steps 10 --optimization
```

Default options are:

```
--steps=50
--optimization=false
```