module ExpInt(
 sf_expint_ei,
 sf_expint_en,
 )
where

import Exp
import Gamma
import Util

----------------------------------------

sf_expint_ei :: Value -> Value
sf_expint_ei z
  | z<0.0  = (0/0)
  | z==0.0 = (-1/0)
  | z<40   = expint_ei__series z
  | otherwise = expint_ei__asymp z

expint_ei__series :: Value -> Value
expint_ei__series z =
  let tterms = ixiter 2 z $ \n t -> t*z/(#)n
      terms = map (\(t,n)->t/(#)n) $ zip tterms [1..]
      res = kahan_sum terms
  in if z<0.5
     then sf_log(z * sf_exp(euler_gamma + res))
     else res + sf_log(z) + euler_gamma

expint_ei__asymp :: Value -> Value
expint_ei__asymp z =
  let terms = tk $ ixiter 1 1.0 $ \n t -> t/z*(#)n
      res = kahan_sum terms
  in res * (sf_exp z) / z
  where tk (a:b:cs) = if (abs a)<(abs b) then [a] else a:(tk$b:cs)

----------------------------------------


sf_expint_en :: Int -> Value -> Value
sf_expint_en n z | z<0 = (0/0) -- TODO: confirm this
                 | z==0 = (1/(#)(n-1)) -- TODO: confirm this
sf_expint_en 0 z = sf_exp(-z) / z
sf_expint_en 1 z = expint_en__1 z
sf_expint_en n z | z<=1.0    = expint_en__series n z
                 | otherwise = expint_en__contfrac n z

expint_en__1 :: Value -> Value
expint_en__1 z =
  let r0 = -euler_gamma - (sf_log z)
      tterms = ixiter 2 (z) $ \k t -> -t*z/(#)k
      terms = map (\(t,k)->t/(#)k) $ zip tterms [1..]
  in kahan_sum (r0:terms)

-- assume n>=2, z<=1
expint_en__series :: Int -> Value -> Value
expint_en__series n z =
  let n' = (#)n
      res = (-(sf_log z) + (sf_digamma n')) * (-z)^(n-1)/(#)(factorial$n-1) + 1/(n'-1)
      terms' = ixiter 2 (-z) (\m t -> -t*z/(#)m)
      terms = map (\(m,t)->(-t)/(#)(m-(n-1))) $ filter ((/=(n-1)).fst) $ zip [1..] terms'
  in kahan_sum (res:terms)

-- assume n>=2, z>1
-- modified Lentz algorithm
expint_en__contfrac :: Int -> Value -> Value
expint_en__contfrac n z =
  let fj = zeta
      cj = fj
      dj = 0
      j = 1
      n' = (#)n
  in lentz j cj dj fj
  where
    zeta = 1e-100
    eps = 5e-16
    nz x = if x==0 then zeta else x
    lentz j cj dj fj =
      let aj = (#) $ if j==1 then 1 else -(j-1)*(n+j-2)
          bj = z + (#)(n + 2*(j-1))
          dj' = nz $ bj + aj*dj
          cj' = nz $ bj + aj/cj
          dji = 1/dj'
          delta = cj'*dji
          fj' = fj*delta
      in if (abs$delta-1)<eps
         then fj' * sf_exp(-z)
         else lentz (j+1) cj' dji fj'

