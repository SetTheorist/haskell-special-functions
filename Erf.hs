{-# Language BangPatterns #-}

module Erf (
    sf_erf,
    sf_erfc,
    erfc_asymp_pos,
    erfc_asymp_pos',
) where

import Exp
import Util

----------------------------------------
-- not perfect, but workable for now

sf_erf :: (Value v) => v -> v
sf_erf z 
  | (re z)<(-1) = -sf_erf(-z)
  | (rabs z)<1  = erf_series z
  | otherwise   = 1 - sf_erfc z

-- infinite loop when (re z)==0
sf_erfc :: (Value v) => v -> v
sf_erfc z 
  | (re z)<(-1) = 2-(sf_erfc (-z))
  | (rabs z)<1  = 1-(sf_erf z)
  | (rabs z)<10 = erfc_cf_pos1 z
  | otherwise   = erfc_asymp_pos z -- TODO: hangs for very large input

erf_series z =
  let z2 = z^2
      rts = ixiter 1 z $ \n t -> (-t)*z2/(#)n
      terms = map (\(n,t)->t/(#)(2*n+1)) $ zip [0..] rts
  in (2/sf_sqrt pi)  * ksum terms

erfc_asymp_pos z =
  let z2 = z^2
      iz2 = 1/2/z2
      terms = ixiter 1 (1/z) $ \n t -> (-t*iz2)*(#)(2*n-1)
      tterms = tk terms
  in (sf_exp (-z2))/(sqrt pi) * ksum tterms
  where tk (a:b:cs) = if (rabs a)<(rabs b) then [a] else a:(tk$b:cs)

-- alternative styling:
erfc_asymp_pos' z =
  let !z2 = z^2
      !iz2 = 1/(2*z2)
      !terms = ixiter 1 (1/z) $ \ !n !t -> (-t*iz2)*(#)(2*n-1)
  in (sf_exp (-z2))/(sqrt pi) * ks terms 0 0
  where ks (t:t':ts) !s !e = kadd t s e $ if rabs(t')>rabs(t) then const else ks (t':ts)

erfc_cf_pos1 z = 
  let z2 = z^2
      as = z:(map fromDouble [1/2,1..])
      bs = 0:cycle [z2,1]
      cf = steeds as bs
  in sf_exp(-z2) / (sqrt pi) * cf

erfc_cf_pos2 z = 
  let z2 = z^2
      as = (2*z):(map (\n->(#)$ -(2*n+1)*(2*n+2)) [0..])
      bs = 0:(map (\n->2*z2+(#)4*n+1) [0..])
      cf = steeds as bs
  in sf_exp(-z2) / (sqrt pi) * cf

