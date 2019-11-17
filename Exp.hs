{-# Language BangPatterns #-}

module Exp (
    sf_exp,
    sf_exp_m1,
    sf_exp_m1vx,
    sf_exp_men,
    sf_exp_menx,
    sf_log,
    sf_log_p1,
) where

import Control.Monad

import Util

----------------------------------------
-- exp(x)
sf_exp :: Value -> Value 
-- exp(-x) = 1/exp(x)
sf_exp x | isInfinite x = if x<0 then 0 else (1.0/0.0)
sf_exp x | x<0 = 1/(sf_exp (-x))
sf_exp x = kahan_sum $ ixiter 1 1.0 $ \n t -> t*x/((#)n)

test__sf_exp :: IO ()
test__sf_exp = do
  cs <- readFile "exp_1.dat"
  let ls = lines cs
  let ws = map words ls
  let ns = map (map read) ws :: [[Value]]
  let re = map (\(x:fx:_)->relerr fx (sf_exp x)) ns
  mapM_ (putStrLn.show) re

----------------------------------------
-- exp(x)-1
sf_exp_m1 :: Value -> Value 
-- exp(-x)-1 = -exp(-x)*(exp(x)-1)
sf_exp_m1 x | isInfinite x = if x<0 then -1 else (1/0)
            | x<0      = -sf_exp x * sf_exp_m1 (-x)
            |otherwise = kahan_sum $ ixiter 2 x $ \n t -> t*x/((#)n)

----------------------------------------
-- (exp(x)-1)/x
sf_exp_m1vx :: Value -> Value
sf_exp_m1vx x
  | isInfinite x = if x<0 then 0 else (1/0)
  | x/=x = x
  | abs(x)>(1/2) = (sf_exp x - 1)/x -- inaccurate for some complex points
  | otherwise =
      let x2 = x^2
      in 2/(2 - x + x2/6/(1
          + x2/(4*(2*3-3)*(2*3-1))/(1
          + x2/(4*(2*4-3)*(2*4-1))/(1
          + x2/(4*(2*5-3)*(2*5-1))/(1
          + x2/(4*(2*6-3)*(2*6-1))/(1
          + x2/(4*(2*7-3)*(2*7-1))/(1
          + x2/(4*(2*8-3)*(2*8-1))/(1
          ))))))));

----------------------------------------

-- Compute tail of series expansion of exponential $(e^z - \sum_(k=0)^(n-1) x^k/k!)$
-- ($n=0, 1, 2, ...$)
sf_exp_men :: Int -> Value -> Value
sf_exp_men n x = (sf_exp_menx n x) * x^n / ((#)$fac n 1)
    where fac 0 !a = a
          fac 1 !a = a
          fac n !a = fac (n-1) (a*n)

-- Compute scaled tail of series expansion of exponential $(e^z - \sum_(k=0)^(n-1) x^k/k!) / (x^n/n!)$
-- ($n=0, 1, 2, ...$)
sf_exp_menx :: Int -> Value -> Value
sf_exp_menx 0 x = sf_exp x
sf_exp_menx 1 x = sf_exp_m1vx x
sf_exp_menx n x
  | isInfinite x = if x>0 then (1/0) else (0) -- TODO: verify
  | x/=x         = x
  | otherwise    = exp_menx__contfrac n x

-- modified Lentz for continued fraction
exp_menx__contfrac n z =
  let !fj = (#)$ n+1
      !cj = fj
      !dj = 0
      !j  = 1
  in lentz j dj cj fj
  where
    !zeta = 1e-150
    !eps = 1e-16
    nz !z = if z==0 then zeta else z
    lentz !j !dj !cj !fj =
      let aj = if (odd j)
               then z*((#)$(j+1)`div`2)
               else -z*((#)$(n+(j`div`2)))
          bj = (#)$n+1+j
          dj' = nz$ bj + aj*dj
          cj' = nz$ bj + aj/cj
          dji = 1/dj'
          deltaj = cj'*dji
          fj' = fj*deltaj
      in if (abs(deltaj-1)<eps)
         then 1/(1-z/fj')
         else lentz (j+1) dji cj' fj'

----------------------------------------

sf_log = log

-- log(1+x)
sf_log_p1 z
  | z/=z = z
  | (abs z)>0.25 = sf_log (1+z)
  | otherwise = log_p1__series z

-- ln(1+x) = 2\sum_{n=0}^\infty (x / x+2)^(2n+1) / 2n+1
log_p1__series z =
  let r = z/(z+2)
      zr2 = r^2
      tterms = iterate (*zr2) (r*zr2)
      terms = map (\(n,t)->t/((#)$2*n+1)) (zip [1..] tterms)
  in 2*(kahan_sum (r:terms))

----------------------------------------

-- $$\ln(1+z) = z/(1+ z/(2+ z/(3+ 4z/(4+ 4z/(5+ 9z/(6+ 9z/(7+ ...)))))))$$
-- NB not bad convergence
ln_1_z_cf z = steeds (z:(ts 1)) [0..]
  where ts n = (n^2*z):(n^2*z):(ts (n+1))

