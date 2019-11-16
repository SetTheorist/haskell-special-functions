{-# Language BangPatterns #-}

--module Quad where
import System.Environment(getArgs)

data Quad = Quad { hi_::Double, lo_::Double }
  deriving (Eq,Show,Ord,Read)

--make a b = qtsum a b Quad
make a b = tsum a b Quad

-- requires |a|>=|b|
qtsum !a !b !k = k s e
  where !s = a+b
        !e = b+(a-s)
        -- !e = b-(s-a)

-- general
tsum !a !b !k = k s e
  where !s = a+b
        !v = s-a
        !e = (a+(v-s))+(b-v)
        -- !e = (a-(s-v))+(b-v)

split !a !k = k ahi alo
  where !t = 134217729*a
        !ahi = t-(t-a)
        !alo = a-ahi

tprod !a !b !k =
  split a $ \ !ahi !alo ->
  split b $ \ !bhi !blo ->
  let !p = a*b
      !e = (((ahi*bhi - p) + ahi*blo) + alo*bhi) + alo*blo
  in k p e

qadd !(Quad !xhi !xlo) !(Quad !yhi !ylo) =
  tsum  xhi yhi     $ \ !hs !he ->
  tsum  xlo ylo     $ \ !ls !le ->
  qtsum hs  (he+ls) $ \ !h  !k  ->
  qtsum h   (le+k)  $ \ !hi !lo ->
  Quad hi lo -- make?

qprod !(Quad !xhi !xlo) !(Quad !yhi !ylo) =
  tprod xhi yhi $ \ p e ->
  qtsum p (e + (xhi*ylo + xlo*yhi)) $ \ hi lo ->
  Quad hi lo -- make?

{--
main :: IO ()
main = do
  args <- getArgs
  let l = length args
  let f = if l==0 then qadd else qadd'
  let s = foldl (\q x->f (make x 0.0) q) (make 0.0 0.0) $ [1/fromInteger n|n<-[1..10000000]]
  putStrLn . show $ s
--}

