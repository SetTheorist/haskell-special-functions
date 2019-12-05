\chapter{Lambert W}

Compute Lambert's W function $W(z) \exp(W(z)) = z$
with optional specification of the branch to return.
(Currently only supports the two real branches b=0, b=-1; and only supports real z).

\begin{code}
{-# Language AllowAmbiguousTypes #-}  -- remove this after cleaning up "Value"
{-# Language BangPatterns #-}
{-# Language FlexibleContexts #-}
{-# Language ScopedTypeVariables #-}
module Lambert where
import Exp
import Util
\end{code}

\section{Lambert}

Positive real branch
\begin{code}
sf_lambert_w :: forall v . (Value v) => (RealKind v) -> (RealKind v)
sf_lambert_w !z 
  | z<(-sf_exp(-1)) = nan
  | z==0      = 0
  | z<0       = lambert__halley z $ -0.1
  | otherwise = lambert__halley z $ sf_log(z / (sf_log_p1 z))
\end{code}

Negative real branch
\begin{code}
sf_lambert_w1 :: forall v . (Value v) => (RealKind v) -> (RealKind v)
sf_lambert_w1 !z 
  | z<(-sf_exp(-1)) || z>0 = nan
  | z==0      = 0
  | z==(-sf_exp(-1)) = -1
  | z<-0.183939 = -- series expansion (for z near -1/e)
      let !p = -sf_sqrt(2*((sf_exp 1)*z + 1))
      in  lambert__halley z $ -1 + p - p^2/3 + p^3*11/72 - 43*p^4/540
  | otherwise = -- asymptotic near 0^-
      let !l1 = sf_log(-z)
          !l2 = sf_log(-sf_log(-z))
      in lambert__halley z $ l1 - l2 + l2/l1 + ((-2+l2)*l2)/(2*l1^2) + ((6-9*l2+2*l2^2)*l2)/(6*l1^3)
\end{code}

Halley iteration
\begin{code}
lambert__halley :: forall v . (Value v, Value(RealKind v)) => (RealKind v) -> (RealKind v) -> (RealKind v)
lambert__halley !z !w = go w
  where 
    go !w =
      let !ew = sf_exp w
          !w' = w - (w*ew - z) / ((w+1)*ew - (w+2)*(w*ew-z)/(2*w+2))
      in if w==w' then w
         else go w'
\end{code}
