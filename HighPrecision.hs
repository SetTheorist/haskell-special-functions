{-# Language BangPatterns #-}

{--
based on
"Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates"
Jonathan Richard Shewchuk (October 1, 1997)
--}

module HighPrecision where

----------------------------------------

-- requires that |a|>=|b|
-- returns non-overlapping expansion x>y such that a+b=x+y
fast_two_sum :: Double -> Double -> (Double -> Double -> a) -> a
fast_two_sum !a !b !k =
  let !x  = a + b
      !bv = x - a
      !y  = b - bv
  in k x y

-- returns non-overlapping expansion x>y such that a+b=x+y
two_sum :: Double -> Double -> (Double -> Double -> a) -> a
two_sum !a !b !k =
  let !x  = a + b
      !bv = x - a
      !av = x - bv
      !br = b - bv
      !ar = a - av
      !y  = ar + br
  in k x y

(^:^) :: Double -> [Double] -> [Double]
(^:^) !0 !xs = xs
(^:^) !x !xs = x:xs

-- requires es is n-term non-overlapping increasing-magnitude expansion
-- returns n+1-term non-overlapping increasing-magnitude expansion hs
-- (that there may be zero terms in input or output)
-- if es is non-adjacent then returned hs is non-adjacent
--   (if round-to-even tiebreaking fp arithmetic is used)
grow_expansion :: [Double] -> Double -> [Double]
grow_expansion !es !b = iter b es
  where iter !q ![]     = q:[]
        iter !q (!e:es) = two_sum q e $ \ !qi !hi ->
          hi:(iter qi es)

-- requires es, fs non-overlapping increasing-magnitude expansions
-- returns n+m-term non-overlapping increasing-magnitude expansion
-- (that there may be zero terms in input or output)
-- if es and fs are non-adjacent then returned hs is non-adjacent
--   (if round-to-even tiebreaking fp arithmetic is used)
expansion_sum :: [Double] -> [Double] -> [Double]
expansion_sum !es !fs = iter [] es fs
  where iter !xs ![]     ![]     = reverse xs
        iter !xs (!h:hs) ![]     = iter (h:xs) hs []
        iter !xs !hs     (!f:fs) =
          let (!hi:hs') = grow_expansion hs f
          in iter (hi:xs) hs' fs

{--
-- TODO: not sure about this...
-- if es,fs are strongly non-overlapping expansions 
-- then result is strong non-overlapping expansion
-- as long as round-to-even rounding is used
fast_expansion_sum :: [Double] -> [Double] -> [Double]
fast_expansion_sum !es !fs =
    let !(g1:g2:gs) = merge es fs
    in fast_two_sum g2 g1 $ \ q2 h1 ->
       iter [h1] q2 gs
  where merge :: [Double] -> [Double] -> [Double]
        merge ![] !bs = bs
        merge !as ![] = as
        --merge !(0:as) !bl = merge as bl
        --merge !al !(0:bs) = merge al bs
        merge !al@(a:as) !bl@(b:bs) = if abs(a)<abs(b) then a:(merge as bl) else b:(merge al bs)
        iter !hs !qim1 ![] = reverse (qim1:hs)
        iter !hs !qim1 !(gi:gs) = two_sum qim1 gi $ \ qi him1 ->
                               iter (him1:hs) qi gs
--}

----------------------------------------

-- hard-coded to split at s=27
-- produces non-overlapping split |alo|<=|ahi|, a=alo+ahi
split :: Double -> (Double -> Double -> a) -> a
split !a !k =
  let !c = 134217729*a
      !abig = c - a
      !ahi = c - abig
      !alo = a - ahi
  in k ahi alo

-- produces non-overlapping expansion such that a*b=x+y, x>y
-- if round-to-even tiebreaking, then x and y are non-adjacent
two_product :: Double -> Double -> (Double -> Double -> a) -> a
two_product !a !b !k =
  split a $ \ !ahi !alo ->
  split b $ \ !bhi !blo ->
  let !x = a*b
      !e1 = x - (ahi*bhi)
      !e2 = e1 - (alo*bhi)
      !e3 = e2 - (ahi*blo)
      !y = (alo*blo) - e3
  in k x y

-- requires non-overlapping, increasing expansion
-- produces non-overlapping, increasing expansion
scale_expansion :: [Double] -> Double -> [Double]
scale_expansion !(e:es) !b = two_product e b $ \ q2 h1 -> h1:(iter q2 es)
  where iter !q ![] = [q]
        iter !q2i2 (!e:es) =
          two_product e b       $ \ !tti  !ti   ->
          two_sum q2i2 ti       $ \ !q2i1 !h2i2 ->
          fast_two_sum tti q2i1 $ \ !q2i  !h2i1 ->
          h2i2:(h2i1:(iter q2i es))

-- produces non-overlapping (non-adjacent if round-to-even)
-- largest component hn approximates h with error < ulp(un)
compress :: [Double] -> [Double]
compress !es =
  let (!f:fs) = reverse es
      (!g:gs) = iterdown f fs
      !hs = iterup g gs
  in reverse hs
  where
    iterdown !q ![] = q^:^[]
    iterdown !q (!f:fs) =
      fast_two_sum q f $ \ !qq' !q' ->
      qq':iterdown q' fs
    iterup q ![] = q^:^[]
    iterup q (!gi:gs) =
      fast_two_sum q gi $ \ !qq' !q' ->
      qq' ^:^ iterup q' gs

----------------------------------------

multiply_expansions :: [Double] -> [Double] -> [Double]
multiply_expansions !es ![] = []
multiply_expansions !es (!f:fs) =
  let xs = scale_expansion es f
  in expansion_sum xs (multiply_expansions es fs)

----------------------------------------

data High = High Int ![Double]
  deriving (Show,Eq,Ord)

make_high n hs = High n $ reverse $ take n $ reverse $ (take n (repeat 0.0))++hs

high_sum n hs = make_high n $ foldl grow_expansion [] hs

(+^) :: High -> Double -> High
(High n es) +^ d = make_high n $ compress $ grow_expansion es d

(-^) :: High -> Double -> High
(High n es) -^ d = make_high n $ compress $ grow_expansion es (-d)

(*^) :: High -> Double -> High
(High n es) *^ d = make_high n $ compress $ scale_expansion es d

-- (/^) :: High -> Double -> High
-- (High n es) /^ d = make_high n $ compress $ ...

showdigs (High n es) = undefined

instance Num High where
  (High n1 h1) + (High n2 h2) = make_high (max n1 n2) $ compress $ expansion_sum h1 h2
  (High n1 h1) - (High n2 h2) = make_high (max n1 n2) $ compress $ expansion_sum h1 (map negate h2)
  (High n1 h1) * (High n2 h2) = make_high (max n1 n2) $ compress $ multiply_expansions h1 h2
  negate (High n1 h1) = make_high n1 $ map negate h1
  abs h@(High n1 h1) = if last h1 < 0 then (negate h) else h
  signum (High n1 h1) = make_high n1 $ [signum $ last h1]
  fromInteger k = High 1 [fromInteger k]

{--
instance Fractional High where
  (/) :: a -> a -> a
  recip :: a -> a
  fromRational :: Rational -> a
--}


