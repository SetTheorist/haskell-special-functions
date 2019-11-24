module Spence (
    sf_spence,
) where

import Exp
import Util


pi2_6 :: (Value v) => v
pi2_6 = pi^2/6

-- Compute Spence's integral.
-- TODO: UNTESTED
-- TODO: verify, especially, complex cases
sf_spence :: (Value v) => v -> v
sf_spence z
  | is_nan z = z
  | (re z)<0   = 0/0
  | z==0  = pi2_6
  | (rabs z)<0.5 =
      (spence__series z) + (pi2_6 - (sf_log z)*(sf_log (1-z)))
  | (rabs z)<1.0 =
      -(spence__series (1-z))
  | (rabs z)<2.5 =
      (spence__series ((z-1)/z)) - (sf_log z)^2/2
  | otherwise = 
      (spence__series (1/(1-z))) - pi2_6 - (sf_log (z-1))^2/2

spence__series z = 
  let zk = iterate (*z) z
      terms = map (\(t,k)-> -t/(#)k^2) (zip zk [1..])
  in ksum terms

