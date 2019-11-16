{-# Language BangPatterns #-}

module Quad where
--import System.Environment(getArgs)

import Data.Char(digitToInt,intToDigit)

data Quad = Quad { hi_::Double, lo_::Double }
  deriving (Eq,Ord,Read)

pi_q :: Quad
pi_q = stoq "3.14159265358979323846264338327950288419716939937510"
eulergamma_q :: Quad
eulergamma_q = stoq "0.57721566490153286060651209008240243104215933593992"
ln2_q :: Quad
ln2_q = stoq "0.69314718055994530941723212145817656807"

rawshow (Quad hi lo) =
  "{"++(show hi)++" "++(show lo)++"}"

--quick and dirty for now
stoq :: String -> Quad
stoq [] = make (0.0/0.0) (0.0/0.0)
stoq ('-':as) = -(stoq as)
stoq ('+':as) = stoq as
stoq as = bef (make 0 0) as
  where bef q [] = q
        bef q ('.':as) = aft q 0 as
        bef q ('e':as) = scale10 (read as) q
        bef q (a:as) =
          let !q' = (q*10) `qdadd` (fromIntegral.digitToInt$a)
          in bef q' as
        aft q n [] = (scale10 n q)
        aft q n ('e':as) = scale10 (n-(read as)) q
        aft q n (a:as) =
          let !q' = (q*10) `qdadd` (fromIntegral.digitToInt$a)
          in aft q' (n+1) as

scale10 :: Int -> Quad -> Quad
scale10 !b !q
  | b>0  = qddivide (scale10 (b-1) q) 10
  | b<0  = qdprod   (scale10 (b+1) q) 10
  | otherwise = q
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
        !q' = scale10 fb q
    in (digs 34 True q')++"e"++(show fb)

instance Show Quad where
  show = fixshow 33

--instance Read Quad where
--  read = stoq

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

-- basic implementation: range-reduce & series
qexp q
  | q<0 = 1/(qexp (-q))
  | otherwise =
    let !nd = q / ln2_q
        !nn = floor.hi_ $ nd
        !r = q - ln2_q * (fromIntegral nn)
        !s = sm r 0.0 1.0 1
    in Quad (scaleFloat nn $ hi_ s) (scaleFloat nn $ lo_ s)
  where sm !r !s !t !n
          | s+t==s = s
          | otherwise = sm r (s+t) (t*r/(fromIntegral n)) (n+1)

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

{--
main :: IO ()
main = do
  args <- getArgs
  let l = length args
  let s = foldl (\q x->(flip qdadd) x q) (make 0.0 0.0) $ [fromInteger (333323*n)|n<-[1..10000000]]
  putStrLn . show $ s
--}

