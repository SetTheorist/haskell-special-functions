module Exp (
    sf_exp,
    sf_expm1,
    sf_exp_m1vx,
    sf_log,
    --sf_log1p,
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
sf_expm1 :: Value -> Value 
-- exp(-x)-1 = -exp(-x)*(exp(x)-1)
sf_expm1 x | isInfinite x = if x<0 then -1 else (1/0)
           | x<0      = -sf_exp x * sf_expm1 (-x)
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
sf_log = log

----------------------------------------

-- $$\ln(1+z) = z/(1+ z/(2+ z/(3+ 4z/(4+ 4z/(5+ 9z/(6+ 9z/(7+ ...)))))))$$
-- NB not bad convergence
ln_1_z_cf z = steeds (z:(ts 1)) [0..]
    where ts n = (n^2*z):(n^2*z):(ts (n+1))

