# Wavelength Effects Simulations

Simulate 100 steps of 10K electrons each, incident uniformly over the central pixel
with an initial electron position randomly sampled from a distribution appropriate
for `FilterBand`.

Run jobs using:
```
nohup src/Poisson data/bands/u-band.cfg > data/bands/u-band.log 2>&1 &
nohup src/Poisson data/bands/g-band.cfg > data/bands/g-band.log 2>&1 &
nohup src/Poisson data/bands/r-band.cfg > data/bands/r-band.log 2>&1 &
nohup src/Poisson data/bands/i-band.cfg > data/bands/i-band.log 2>&1 &
nohup src/Poisson data/bands/z-band.cfg > data/bands/z-band.log 2>&1 &
nohup src/Poisson data/bands/y-band.cfg > data/bands/y-band.log 2>&1 &
```

Make plots using:
```
export POISSON_BANDS=/data/desc/poisson/bands
cd data/bands
./make_plots.py
```
