from __future__ import division
import numpy as np

# 1.65: 95% (colored if (3,0)
# 1.25: 90% (colored if (2,0)
def intervals(u,d,z=1.65):
  a = 1. + u
  b = 1. + d
  mu_a = a/(a+b)
  mu_b = b/(a+b)
  std_err = z*np.sqrt( (a*b)/( (a+b)**2*(a+b+1.) ) )
  return ( mu_a, mu_b, std_err )

# subtract std_err from mean. If over 50%, add that much color.

# color by most highly enriched, tint by confidence
# blend based on three quads


# what if dimension is continuous, like time? Multiple classes?
# TODO:
# ------------------------------
# compute all pairs boolean, count number of points in each quad by class
