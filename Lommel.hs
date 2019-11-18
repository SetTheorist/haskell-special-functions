module Lommel (
  sf_lommel_s,
  sf_lommel_s2,
) where

import Util

--TODO: These are completely untested!

-- Compute the Lommel function $s_(\mu,\nu)(z)$
sf_lommel_s mu nu z = 
  -- series (breaks down when mu+/-nu is odd integer...)
  let terms = ixiter 1 1.0 $ \ k t -> -t*z^2 / ((mu+((#)$2*k+1))^2 - nu^2)
      res = kahan_sum terms
  in res * z**(mu+1) / ((mu+1)^2 - nu^2)


-- Compute the (asymptotic) second Lommel function $S_(\mu,\nu)(z) ~ \sum ...$
sf_lommel_s2 mu nu z =
  -- asymptotic series (breaks down when mu+/-nu is odd integer...)
  let tterms = ixiter 1 1.0 $ \ k t -> -t*((mu-((#)$2*k+1))^2 - nu^2) / z^2
      terms = tk tterms
      res = kahan_sum terms
  in res
  where tk (a:b:cs) = if abs(a)<abs(b) then [a] else a:(tk$b:cs)
