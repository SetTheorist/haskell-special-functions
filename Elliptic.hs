module Elliptic where

import AGM
import Util

----------------------------------------

-- Compute the complete elliptic integral of the first kind $K(k)$
sf_elliptic_k :: Value -> Value
sf_elliptic_k k =
  let an = sf_agm' 1.0 (sf_sqrt $ 1.0-k^2)
  in pi/(2*an)

----------------------------------------

-- Compute the symmetric elliptic integral of the first kind $R_F(x,y,z)$
-- x,y,z>0
sf_elliptic_rf :: Double -> Double -> Double -> Double
sf_elliptic_rf x y z = 1/(sf_sqrt $ iter x y z)
  where
    iter x y z =
      let lam = (sf_sqrt $ x*y) + (sf_sqrt $ y*z) + (sf_sqrt $ z*x)
          mu = (x+y+z)/3
          eps = foldl1 max $ map (\a->abs(1-a/mu)) [x,y,z]
          x' = (x+lam)/4
          y' = (y+lam)/4
          z' = (z+lam)/4
      in if (eps<1e-16) || ([x,y,z]==[x',y',z'])
         then x
         else iter x' y' z'

