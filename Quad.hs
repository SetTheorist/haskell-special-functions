{-# Language BangPatterns #-}

--module Quad where
import System.Environment(getArgs)

import Data.Char(intToDigit)

data Quad = Quad { hi_::Double, lo_::Double }
  deriving (Eq,Ord,Read)

rawshow (Quad hi lo) =
  "{"++(show hi)++" "++(show lo)++"}"


scale :: Int -> Quad -> Quad
scale !b !q
  | b==0 = q
  | b>0  = qddivide (scale (b-1) q) 10
  | b<0  = qdprod   (scale (b+1) q) 10
digs !n !p !q
  | n==0 = []
  | otherwise =
      let !d = max 0 $ floor (hi_ q)
          !s = (intToDigit d):(if p then "." else "")
          !q' = qdadd q (fromIntegral (-d)) `qdprod` 10
      in s++(digs (n-1) False q')
fixshow n q@(Quad hi lo)
  | q<(Quad 0 0) = '-':(fixshow n (qnegate q))
  | otherwise = 
    let !b  = logBase 10 hi
        !fb = floor b
        !q' = scale fb q
    in (digs 34 True q')++"e"++(show fb)

instance Show Quad where
  show = fixshow 33

--make a b = qtsum a b Quad
make :: Double -> Double -> Quad
make !a !b = ddsum a b Quad

parts :: Quad -> (Double -> Double -> a) -> a
parts !(Quad !hi !lo) k = k hi lo

-- requires |a|>=|b|
qtsum :: Double -> Double -> (Double -> Double -> a) -> a
qtsum !a !b !k = k s e
  where !s = a+b
        !e = b+(a-s)
        -- !e = b-(s-a)

-- general
ddsum :: Double -> Double -> (Double -> Double -> a) -> a
ddsum !a !b !k = k s e
  where !s = a+b
        !v = s-a
        !e = (a+(v-s))+(b-v)
        -- !e = (a-(s-v))+(b-v)

split :: Double -> (Double -> Double -> a) -> a
split !a !k = k ahi alo
  where !t = 134217729*a
        !ahi = t - (t - a)
        !alo = a - ahi

ddprod :: Double -> Double -> (Double -> Double -> a) -> a
ddprod !a !b !k =
  split a $ \ !ahi !alo ->
  split b $ \ !bhi !blo ->
  let !p = a*b
      !e = (((ahi*bhi - p) + ahi*blo) + alo*bhi) + alo*blo
  in k p e

qdadd :: Quad -> Double -> Quad
qdadd !(Quad !xhi !xlo) y = 
  ddsum y   xhi         $ \ shi slo ->
  qtsum shi (slo + xlo) $ \ hhi hlo ->
  qtsum hhi hlo         $ \ hi  lo -> 
  Quad hi lo
dqadd :: Double -> Quad -> Quad
dqadd = flip qdadd

qqadd :: Quad -> Quad -> Quad
qqadd !(Quad !xhi !xlo) !(Quad !yhi !ylo) =
  ddsum  xhi yhi    $ \ !hs !he ->
  ddsum  xlo ylo    $ \ !ls !le ->
  qtsum hs  (he+ls) $ \ !h  !k  ->
  qtsum h   (le+k)  $ \ !hi !lo ->
  Quad hi lo

qdprod :: Quad -> Double -> Quad
qdprod !(Quad !xhi !xlo) !y =
  ddprod xhi y           $ \ !thi !tlo ->
  qtsum  thi (tlo+y*xlo) $ \ !hi  !lo ->
  Quad hi lo

qqprod :: Quad -> Quad -> Quad
qqprod !(Quad !xhi !xlo) !(Quad !yhi !ylo) =
  ddprod xhi yhi $ \ p e ->
  qtsum p (e + (xhi*ylo + xlo*yhi)) $ \ hi lo ->
  Quad hi lo

qnegate :: Quad -> Quad
qnegate !(Quad !hi !lo) = Quad (negate hi) (negate lo)

qqdivide :: Quad -> Quad -> Quad
qqdivide !(Quad !xhi !xlo) !(Quad !yhi !ylo) =
  let !cc = xhi / yhi
  in ddprod cc yhi $ \ !uu !u ->
     let !c = ((((xhi-uu)-u)+xlo)-cc*ylo)/yhi
     in qtsum cc c $ \ !hi !lo ->
        Quad hi lo

dqdivide :: Double -> Quad -> Quad
dqdivide !x !(Quad !yhi !ylo) =
  let !cc = x / yhi
  in ddprod cc yhi $ \ !uu !u ->
     let !c = ((((x-uu)-u))-cc*ylo)/yhi
     in qtsum cc c $ \ !hi !lo ->
        Quad hi lo

qddivide :: Quad -> Double -> Quad
qddivide !(Quad !xhi !xlo) !y = 
  let !xdy = xhi / y
  in ddprod xdy y $ \ !uu !u ->
     let !c = (((xhi-uu)-u)+xlo)/y
     in qtsum xdy c $ \ !hi !lo ->
        Quad hi lo

instance Num Quad where
  (+) = qqadd
  (-) a b = a + (-b)
  (*) = qqprod
  negate = qnegate
  abs q = if q<0 then (qnegate q) else q
  signum (Quad !hi _) = make (signum hi) 0
  -- TODO: fix for lower-order bits
  fromInteger i = (make (fromIntegral i) 0)

instance Fractional Quad where
  (/) = qqdivide
  --recip :: a -> a
  -- TODO: fix for lower-order bits
  fromRational r = make (fromRational r) 0
       

main :: IO ()
main = do
  args <- getArgs
  let l = length args
  let s = foldl (\q x->(flip qdadd) x q) (make 0.0 0.0) $ [fromInteger (333323*n)|n<-[1..10000000]]
  putStrLn . show $ s

