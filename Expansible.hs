{-# Language BangPatterns #-}

import Quad

class Doublish e where
  (+.) :: e -> Double -> e
  (-.) :: e -> Double -> e
  trunc :: e -> Double
  embed :: Double -> e
  magic :: e -> Int

{--
data Expanded a = Expanded { hi_::a, lo_::Double }
widen :: a -> Expanded a
widen !x = Expanded x 0
narrow :: Expanded a -> a
narrow !x = hi_ x
--}

instance Doublish Double where
  (+.) = (+)
  (-.) = (-)
  trunc = id
  embed = id
  magic _ = 27

instance Doublish Quad where
  q +. d = qdadd q d
  q -. d = qdadd q (-d)
  trunc = hi_
  embed x = (Quad x 0)
  magic _ = 54

add_d :: (Num a, Doublish a) => a -> Double -> (a -> Double -> t) -> t
add_d !a !b !k = k s e
  where !s = a +. b
        !v = s - a
        !e = trunc $ (a+(v-s))+((embed b)-v)

sub_d :: (Num a, Doublish a) => a -> Double -> (a -> Double -> t) -> t
sub_d !a !b !k = add_d a (-b) k

splitt :: (Num a, Doublish a) => a -> (a -> a -> t) -> t
splitt a k = k ahi alo
  where !t = (2^(magic a) + 1)*a
        !ahi = t - (t - a)
        !alo = a - ahi

{--
split :: Double -> (Double -> Double -> a) -> a
split !a !k = k ahi alo
  where !t = 134217729*a
        !ahi = t - (t - a)
        !alo = a - ahi
--}

