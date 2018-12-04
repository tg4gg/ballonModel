# Balloon Model
MBSE MitX Week 2 project - Balloon modelling

This script models the movement of a balloon filled with helium from the bottom
to the ceiling of a building, it uses the forward Euler approximation method.
If we set the flag `refinedModel` it'll recalculate values for air pressure
in each step and will consider a meandering motion path for the balloon.

## What are the dependencies?
It uses numpy and matplotlib, if you're not familiar with these packages 
you can download Anaconda which is available for MacOS, Linux and Windows.

## Whats the output?
It prints the parameters and the result including the error compared to
the experimental time.
```
Lm=27.599, Dm=0.275, m_b_kg=0.0029, V=0.0109, A=0.0594
vt=2.241, m_tot=0.0103, alpha_C=0.0103, beta_C=0.0181, gamma_C=-0.0371

RESULTS: t=21.226, ddz=-0.000, dz=1.317, z=27.599, ind=21227
ERROR: -0.12% Model time=21.226 Experimental time=21.2
```

It'll also print graphs with the data:
![alt balloon model](https://github.com/tg4gg/ballonModel/blob/master/example.png)
Zooming in into the first 2.5 seconds:
![alt ballon model zoomed](https://github.com/tg4gg/ballonModel/blob/master/example-zoom.png)
