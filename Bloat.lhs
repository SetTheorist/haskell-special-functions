\chapter{Big Float}

\section{Preamble}

\begin{code}
{-# Language BangPatterns #-}
{-# Language MagicHash #-}
module Bloat where

import GHC.Exts(Int(I#))
import GHC.Integer.Logarithms(integerLog2#)

import Data.Bits(shift)
import Data.Char(intToDigit)
import Data.Ratio

infix 2 ?
(?) !b x y = if b then x else y
{-# INLINE (?) #-}
\end{code}

\section{Datatypes}

%\subsection{Sign}
\begin{code}
data Sign = Positive | Negative
  deriving (Show,Eq)

mulsign !Positive !Negative = Negative
mulsign !Negative !Positive = Negative
mulsign !_        !_        = Positive

sign2int !Positive = 1
sign2int !Negative = -1
\end{code}

%\subsection{Bloat}
\begin{code}
data Bloat = Bloat {
    bdig :: !Int,
    bsgn :: !Sign,
    bman :: !Integer,
    bexp :: !Int
  }
  deriving (Show,Eq)
\end{code}

\section{Utility functions}

\begin{code}
ilog2 !n = 1 + (I#(integerLog2# n))

binshow !0 = "0"
binshow !n = go n []
  where
    go !0 !acc = acc
    go !n !acc =
      let (q,r) = n`quotRem`2
      in go q ((r==0 ? '0' $ '1'):acc)

-- helper func, assumes 1<=r<10
rshow !_ !0 !acc = reverse acc
rshow !0 !_ !acc = reverse acc
rshow !d !r !acc = rshow d' r' acc'
  where
    !d' = d-1
    !k  = max 0 $ floor r
    !r' = 10 * (r - (fromIntegral k))
    !c  = intToDigit k
    !acc' = c:acc

-- assumes r>=0
r10 !r !e
  | r==0      = (0,0)
  | r>=10     = r10 (r/10) (e+1)
  | r<1       = r10 (r*10) (e-1)
  | otherwise = (r,e)

bshow b | b_isposinf b = "B+Inf"
bshow b | b_isneginf b = "B-Inf"
bshow b | b_isnan b = "BNaN"
bshow b | b_iszero b = "B0"
bshow !(Bloat d s m e) =
  let !r = (toRational m) * 2^^e
      !(r',e') = r10 r 0
      !digs = 1 + (floor $ (fromIntegral d)*0.30102999566398114)
      !rs = rshow digs r' []
  in ((head rs):'.':(tail rs)) ++ ('e':(show e'))
\end{code}

\begin{code}
b_isinf    !b = (bman b)==0 && (bexp b)==1
b_isposinf !b = (b_isinf b) && (bsgn b)==Positive
b_isneginf !b = (b_isinf b) && (bsgn b)==Negative
b_isnan    !b = (bman b)==0 && (bexp b)==(-1)
b_iszero   !b = (bman b)==0 && (bexp b)==0
\end{code}


%\subsection{Constructors}
\begin{code}
bnan    d = Bloat d Positive 0 (-1)
bposinf d = Bloat d Positive 0 1
bneginf d = Bloat d Negative 0 1
bzero   d = Bloat d Positive 0 0

i2b !d !0 = Bloat d Positive 0 0
i2b !d' !n = prec d' $! Bloat d s m e
  where
    !d = ilog2 m
    !s = n>0 ? Positive $ Negative
    !m = abs n
    !e = 0

-- truncates / extends mantissa
-- TODO: do correct rounding here!
prec !d' !(Bloat d s 0 e) = Bloat d' s 0 e
prec !d' !(Bloat d s m e) = Bloat d' s m' e'
  where
    !k  = d' - d
    !m' = m `shift` k
    !e' = e - k

-- extends mantissa by k digits
bshift !k !(Bloat d s m e) = Bloat d' s' m' e'
  where
    !d' = d + k
    !s' = s
    !m' = m `shift` k
    !e' = e - k
\end{code}


\section{Basic Mathematical Operations}

\subsection{Addition, subtraction, multiplication}
\begin{code}
-- any sign, doesn't handle special values
bmul !d' !(Bloat d1 s1 m1 e1) !(Bloat d2 s2 m2 e2) = prec d' $! Bloat d s m e
  where
    !d = ilog2 m
    !s = mulsign s1 s2
    !m = m1 * m2
    !e = e1 + e2

-- assumes identical sign, doesn't handle special values
badd_ss !d' !(Bloat d1 s1 m1 e1) !(Bloat d2 s2 m2 e2) = prec d' $! Bloat d s m e
  where
    !e = min e1 e2
    !m1' = m1 `shift` (e1 - e)
    !m2' = m2 `shift` (e2 - e)
    !d = ilog2 m
    !s = s1
    !m = m1' + m2'

-- assumes opposite sign, doesn't handle special values
badd_ds !d' !(Bloat d1 s1 m1 e1) !(Bloat d2 s2 m2 e2) = prec d' $! (m==0 ? (bzero d') $ (Bloat d s m e))
  where
    !e = min e1 e2
    !m1' = m1 `shift` (e1 - e)
    !m2' = m2 `shift` (e2 - e)
    !d = ilog2 m
    !s = m1'>=m2' ? s1 $ s2
    !m = m1'>=m2' ? (m1' - m2') $ (m2' - m1')

badd !d' b1 b2
  | b_iszero b1 = prec d' b2
  | b_iszero b2 = prec d' b1
  | b_isnan  b1 = prec d' b1
  | b_isnan  b2 = prec d' b2
  | b_isinf  b1 && b_isinf b2 = (bsgn b1==bsgn b2) ? prec d' b1 $ bnan d'
  | b_isinf  b1 = prec d' b1
  | b_isinf  b2 = prec d' b2
  | bsgn b1 == bsgn b2 = badd_ss d' b1 b2
  | otherwise          = badd_ds d' b1 b2

babs  !(Bloat d _ m e) = Bloat d Positive m e

bsign !(Bloat d Positive 0 0) = bzero d
bsign !(Bloat d Positive _ _) = i2b d 1
bsign !(Bloat d Negative _ _) = i2b d (-1)

bnegate !(Bloat d s m e) = Bloat d (s==Positive?Negative$Positive) m e
\end{code}


\begin{code}
-- TODO: handle special cases!
instance Num Bloat where
  b1 + b2 = badd (max (bdig b1) (bdig b2)) b1 b2
  b1 - b2 = b1 + (negate b2)
  b1 * b2 = bmul (max (bdig b1) (bdig b2)) b1 b2
  negate = bnegate
  abs = babs
  signum = bsign

  -- Maybe have a "default digits" setting for this?
  -- Auto-expand in operations (though what about division??!)
  fromInteger 0 = bzero 1
  fromInteger n | n>0 = i2b (ilog2 n) n
  fromInteger n | n<0 = negate.fromInteger.abs $ n
\end{code}

\subsection{Division}

\begin{code}
-- no special cases handled
-- assumes d1 >> d2
bdiv'' d' !(Bloat d1 s1 m1 e1) !(Bloat d2 s2 m2 e2) = prec d' $ Bloat d s m e
  where !d = ilog2 m
        !s = mulsign s1 s2
        !m = m1 `div` m2
        !e = e1 - e2

-- no special cases handled
bdiv' !d' !b1 !b2 = prec d' $ bdiv'' d' (prec (2*d'+4) b1) (prec d' b2)

-- TODO: handle special cases!
bdiv !d' !b1 !b2 = bdiv' d' b1 b2
\end{code}


\begin{code}
instance Fractional Bloat where
  b1 / b2 = bdiv (max (bdig b1) (bdig b2)) b1 b2
  recip b2 = (i2b 1 1) / b2
  fromRational r = (fromInteger.numerator $ r) / (fromInteger.denominator $ r)
{--
main :: IO ()
main = do
  putStrLn $ bshow $ (i2b 300 1) / 10
--}
\end{code}




