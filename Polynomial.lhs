\chapter{Polynomial}

A basic data class for easier manipulation of polynomials.
This is mainly intended for internal use.

\section{Preamble}

\begin{code}
{-# Language BangPatterns #-}
-- {-# Language InstanceSigs #-}
module Polynomial where
\end{code}

\begin{code}
data Poly a = Poly [a]

pzero :: Poly a
pzero = Poly []

px :: Num a => Poly a
px = Poly [0,1]

coeffs (Poly as) = as

instance (Eq a, Num a, Show a) => Show (Poly a) where
  show (Poly [])     = "0"
  show (Poly (a:as)) = (show a) ++ (foldl (++) "" (zipWith st [1..] as))
    where st n 0 = ""
          st 1 1 = " + x"
          st 1 a = " + "++(show a)++"x"
          st n 1 = " + x^"++(show n)
          st n a = " + "++(show a)++"x^"++(show n)

instance Functor Poly where
  fmap f (Poly as) = Poly $ fmap f as

normalize :: (Eq a, Num a) => Poly a -> Poly a
normalize (Poly as) = Poly $ reverse . dropWhile (==0) . reverse $ as
\end{code}

\begin{code}
zipWithDefault :: a -> b -> (a -> b -> c) -> [a] -> [b] -> [c]
zipWithDefault da db f = z
  where
    z []     []     = []
    z []     (b:bs) = (f da b ):(z [] bs)
    z (a:as) []     = (f a  db):(z as [])
    z (a:as) (b:bs) = (f a  b ):(z as bs)
\end{code}


Scalar-polynomial operations
\begin{code}
(+.) :: (Num a) => a -> [a] -> [a]
(+.) sa pb = fmap (sa+) pb

(.+) :: (Num a) => [a] -> a -> [a]
(.+) pa sb = fmap (+sb) pa

(-.) :: (Num a) => a -> [a] -> [a]
(-.) sa pb = fmap (sa-) pb

(.-) :: (Num a) => [a] -> a -> [a]
(.-) pa sb = fmap (subtract sb) pa

(*.) :: (Num a) => a -> [a] -> [a]
(*.) sa pb = fmap (sa*) pb

(.*) :: (Num a) => [a] -> a -> [a]
(.*) pa sb = fmap (*sb) pa

(/.) :: (Fractional a) => a -> [a] -> [a]
(/.) sa pb = fmap (sa/) pb

(./) :: (Fractional a) => [a] -> a -> [a]
(./) pa sb = fmap (/sb) pa
\end{code}

Polynomial-polynomial operations
\begin{code}
(.+.) :: (Num a) => [a] -> [a] -> [a]
(.+.) = zipWithDefault 0 0 (+)

(.-.) :: (Num a) => [a] -> [a] -> [a]
(.-.) = zipWithDefault 0 0 (-)

(.*.) :: (Num a) => [a] -> [a] -> [a]
_      .*. []     = []
[]     .*. _      = []
(a:[]) .*. (b:bs) = (a*b):          (a*.bs)
(a:as) .*. (b:[]) = (a*b):(as.*b)
(a:as) .*. (b:bs) = (a*b):(as.*b).+.(a*.bs).+.(0:(as.*.bs))
\end{code}

\begin{code}
instance (Eq a, Num a) => Num (Poly a) where
  -- (+) :: Poly a -> Poly a -> Poly a
  (Poly pa) + (Poly pb) = normalize $ Poly $ pa .+. pb
  (Poly pa) - (Poly pb) = normalize $ Poly $ pa .-. pb
  (Poly pa) * (Poly pb) = normalize $ Poly $ pa .*. pb

  negate = fmap negate

  -- TODO: fix this
  abs = fmap abs
  signum = normalize . (fmap signum)

  fromInteger 0 = Poly []
  fromInteger i = Poly [fromInteger i]

(Poly as) +> sb = Poly $ zipWith (+) as [sb]
sa <+ (Poly bs) = Poly $ zipWith (+) [sa] bs

(Poly as) *> sb = Poly $ as .* sb
sa <* (Poly bs) = Poly $ sa *. bs

(Poly as) /> sb = Poly $ as ./ sb
\end{code}



Divides polynomials.
\begin{code}
divide_ :: (Fractional a) => [a] -> [a] -> ([a], [a])
divide_ pn pd = undefined
\end{code}


\begin{code}
derivative (Poly []) = Poly $ []
derivative (Poly as) = Poly $ (zipWith (*) (tail as) [1..])

integral   (Poly []) = Poly []
integral   (Poly as) = Poly $ 0:(zipWith (/) as [1..])
\end{code}


