module Bessel where

import Gamma
import Util

----------------------------------------
-- Bessel J(nu,x)
bessel_j_series :: Value -> Value -> Value
bessel_j_series nu z = 
  let z2 = -(z/2)^2
      terms = ixiter 1 1 $ \n t -> t*z2/((#) n)/(nu+(#) n)
      res = kahan_sum terms
  in res * (z/2)**nu / sf_gamma (1+nu)
