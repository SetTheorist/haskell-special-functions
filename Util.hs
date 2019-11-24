{-# Language BangPatterns #-}
{-# Language FlexibleContexts #-}
{-# Language FlexibleInstances #-}
{-# Language TypeFamilies #-}
-- {-# Language UndecidableInstances #-}

module Util where

import Data.Complex

--------------------------------------------------------------------------------

type CDouble = Complex Double

class (Eq t, Floating t, Fractional t, Num t,
       Enum (RealKind t), Eq (RealKind t), Floating (RealKind t), Fractional (RealKind t), Num (RealKind t), Ord (RealKind t),
       Eq (ComplexKind t), Floating (ComplexKind t), Fractional (ComplexKind t), Num (ComplexKind t)
      ) => Value t where
  type RealKind t :: *
  type ComplexKind t :: *
  re :: t -> (RealKind t)
  im :: t -> (RealKind t)
  rabs :: t -> (RealKind t)
  is_inf :: t -> Bool
  is_nan :: t -> Bool
  fromDouble :: Double -> t
  fromReal :: (RealKind t) -> t
  toComplex :: t -> (ComplexKind t)

instance Value Double where
  type RealKind Double = Double
  type ComplexKind Double = CDouble
  re = id
  im = const 0
  rabs = abs
  is_inf = isInfinite
  is_nan = isNaN
  fromDouble = id
  fromReal = id
  toComplex x = x :+ 0

instance Value CDouble where
  type RealKind CDouble = Double
  type ComplexKind CDouble = CDouble
  re = realPart
  im = imagPart
  rabs = realPart.abs
  is_inf z = (is_inf.re$z) || (is_inf.im$z)
  is_nan z = (is_nan.re$z) || (is_nan.im$z)
  fromDouble x = x :+ 0
  fromReal x = x :+ 0
  toComplex = id

-- TODO: add quad versions also

--------------------------------------------------------------------------------

(#) :: (Integral a, Num b) => a -> b
(#) = fromIntegral

relerr :: (Value v) => v -> v -> (RealKind v)
relerr !ex !ap = re $! logBase 10 (abs ((ex-ap)/ex))

ksum :: (Value v) => [v] -> v
ksum terms = k 0 0 terms
  where
    k !sum !err [] = sum
    k !sum !err (t:terms) =
      let !y = t - err
          !x = sum + y
          !err' = (x - sum) - y
      in if x==sum
         then sum
         else (k x err' terms)


-- kadd value oldsum olderr ---> newsum newerr
kadd :: (Value v) => v -> v -> v -> (v -> v -> a) -> a
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
steeds :: (Value v) => [v] -> [v] -> v
steeds (a1:as) (b0:b1:bs) =
    let !c0 = b0
        !d1 = 1/b1
        !delc1 = a1*d1
        !c1 = c0 + delc1
    in recur c1 delc1 d1 as bs
    --in c0:c1:(recur c1 delc1 d1 as bs)
    where recur !cn_1 !delcn_1 !dn_1 !(an:as) !(bn:bs) = 
            let !dn = 1/(dn_1*an+bn)
                !delcn = (bn*dn - 1)*delcn_1
                !cn = cn_1 + delcn
            in if (cn==cn_1) || is_nan cn then cn else (recur cn delcn dn as bs)

--------------------------------------------------------------------------------

sf_sqrt :: (Value v) => v -> v
sf_sqrt = sqrt

--------------------------------------------------------------------------------

