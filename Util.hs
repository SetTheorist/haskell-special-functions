module Util where

type Value = Double

(#) :: (Integral a, Num b) => a -> b
(#) = fromIntegral

sf_sqrt = sqrt

relerr :: Value -> Value -> Value
relerr ex ap = logBase 10 (abs ((ex-ap)/ex))

kahan_sum :: [Value] -> Value
kahan_sum terms = k 0.0 0.0 terms
  where
    k sum err [] = sum
    k sum err (t:terms) =
      let y = t - err
          x = sum + y
          err' = (x - sum) - y
      in if x==sum
         then sum
         else (k x err' terms)


-- kadd value oldsum olderr ---> newsum newerr
kadd :: Value -> Value -> Value -> (Value -> Value -> a) -> a
kadd t s e k =
  let y = t - e
      s' = s + y
      e' = (s' - s) - y
  in k s' e'

ixiter :: (Enum ix) => ix -> a -> (ix->a->a) -> [a]
ixiter i x f = x:(ixiter (succ i) (f i x) f)

-- Steed's algorithm for evaluation of continued fraction
-- C = b0 + a1/(b1+ a2/(b2+ ...))
-- Cn = ...an/bn = An/Bn
--steeds :: [Value] -> [Value] -> [Value]
steeds :: [Value] -> [Value] -> Value
steeds (a1:as) (b0:b1:bs) =
    let c0 = b0
        d1 = 1/b1
        delc1 = a1*d1
        c1 = c0 + delc1
    in recur c1 delc1 d1 as bs
    --in c0:c1:(recur c1 delc1 d1 as bs)
    where recur cn_1 delcn_1 dn_1 (an:as) (bn:bs) = 
            let dn = 1/(dn_1*an+bn)
                delcn = (bn*dn - 1)*delcn_1
                cn = cn_1 + delcn
            in if (cn==cn_1) then cn else (recur cn delcn dn as bs)

