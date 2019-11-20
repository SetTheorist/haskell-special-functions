module PowerSeries (
    PowerSeries,
    fromList,
    toList,
    Num,
) where

----------------------------------------
-- data definition

infixr 4 :..
data PowerSeries a = Zero | a:..(PowerSeries a)
    deriving (Eq,Show,Ord)

fromList [] = Zero
fromList (a:as) = a:..(fromList as)

toList Zero = []
toList (a:..as) = a:(toList as)

----------------------------------------
-- scalar operations

Zero +. s = s:..Zero
(a:..pa) +. s = (a+s):..pa

Zero -. s = (-s):..Zero
(a:..pa) -. s = (a-s):..pa

Zero *. s = Zero
(a:..pa) *. s = (a*s):..(pa*.s)

Zero /. s = Zero
(a:..pa) /. s = (a/s):..(pa/.s)

----------------------------------------
-- series operations

Zero .+. q = q
p    .+. Zero = p
(a:..pa) .+. (b:..pb) = (a+b):..(pa.+.pb)

Zero .-. q = negate q
p    .-. Zero = p
(a:..pa) .-. (b:..pb) = (a-b):..(pa.-.pb)

Zero .*. q     = Zero
p    .*. Zero  = Zero
(a:..pa) .*. (b:..pb) = (a*b):..((pa*.b).+.(pb*.a).+.(0:..(pa.*.pb)))

----------------------------------------

instance Functor PowerSeries where
    fmap f Zero = Zero
    fmap f (a:..pa) = (f a):..(fmap f pa)

----------------------------------------

instance (Num a) => Num (PowerSeries a) where
  pa + pb = pa .+. pb
  pa - pb = pa .-. pb
  pa * pb = pa .*. pb
  negate = fmap negate
  abs = fmap abs
  signum = fmap signum
  fromInteger i = (fromInteger i):..Zero


