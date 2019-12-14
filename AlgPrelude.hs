{-# Language BangPatterns #-}
{-# Language MagicHash #-}
{-# Language NoImplicitPrelude #-}
{-# Language RebindableSyntax #-}

module AlgPrelude (
  ) where

-- import qualified Prelude as P
import Prelude(Bool,Double,Rational,(==),Ord,(<),error,encodeFloat)
import qualified Prelude as P(fromIntegral,(++))

import GHC.Integer(
  Integer,
  plusInteger, minusInteger, timesInteger,
  )
import GHC.Prim(
  Double#, negateDouble#,
  (+##),(-##),(*##),(/##),
  (<##),(<=##),(>##),(>=##),
  (+#), (-#), (*#), negateInt#,
  )
import GHC.Types(
  Double(D#),
  isTrue#,
  Int(I#),
  )

----------------------------------------

infixl 6 +
infixl 6 -
infixl 7 *
infixl 7 /
infixr 8 ^
infixr 8 ^^
infixr 8 *^
infixr 8 *^^

----------------------------------------

-- TODO: do this the smart way!
-- assumes associativity of operation
power :: (Integral i) => (t -> t -> t) -> t -> (t -> t) -> t -> i -> t
power o i n b e
  | e < 0     = power o i n (n b) (negate e)
  | e == 0    = i
  | e == e    = b `o` (power o i n b (e-1))

----------------------------------------

class AddSemigroup t where
  (+) :: t -> t -> t

class AddSemigroup t => AddMonoid t where
  zero :: t
  (*^) :: (Integral i) => t -> i -> t
  (*^) = power (+) zero (error "Negative multiple")

class AddMonoid t => AddGroup t where
  (-) :: t -> t -> t
  negate :: t -> t
  negate = (zero -)
  (*^^) :: (Integral i) => t -> i -> t
  (*^^) = power (+) zero negate

----------------------------------------

class MulSemigroup t where
  (*) :: t -> t -> t

class MulSemigroup t => MulMonoid t where
  one :: t
  (^) :: (Integral i) => t -> i -> t
  (^) = power (*) one (error "Negative power")

class MulMonoid t => MulGroup t where
  (/) :: t -> t -> t
  recip :: t -> t
  recip = (one /)
  (^^) :: (Integral i) => t -> i -> t
  (^^) = power (*) one (error "Negative power")
  -- (^^) = power (*) one recip

----------------------------------------

class (AddGroup t, MulMonoid t) => Ring t where
  -- zembed :: Integral i => i -> t

class (Ring t, MulGroup t) => Field t where
  iszero :: t -> Bool
  qembed :: Rational -> t

----------------------------------------

class (Ord t, Ring t) => Integral t where
  fromInteger :: Integer -> t
  toInteger :: t -> Integer

class (Ord t, Field t) => FromRational t where
  fromRational :: Rational -> t

----------------------------------------
----------------------------------------

instance AddSemigroup Int where (+) = plusInt
instance AddMonoid    Int where zero = 0
instance AddGroup     Int where (-) = minusInt
instance MulSemigroup Int where (*) = timesInt
instance MulMonoid    Int where one = 1
instance Ring         Int -- where zembed i = P.fromIntegral i
instance Integral     Int where
  fromInteger x = P.fromIntegral x
  toInteger x = P.fromIntegral x

plusInt   (I# x) (I# y) = I# (x +# y)
minusInt  (I# x) (I# y) = I# (x -# y)
negateInt (I# x)        = I# (negateInt# x)
timesInt  (I# x) (I# y) = I# (x *# y)

{--
----------------------------------------
----------------------------------------

instance AddSemigroup Integer where (+) = plusInteger
instance AddMonoid    Integer where zero = 0
instance AddGroup     Integer where (-) = minusInteger
instance MulSemigroup Integer where (*) = timesInteger
instance MulMonoid    Integer where one = 1
instance Ring Integer where zembed i = i
instance Integral Integer where fromInteger x = x

----------------------------------------
----------------------------------------

instance Integral Double where
  fromInteger x = let x = x in x

instance FromRational Double where
  fromRational x = let x = x in x

instance AddSemigroup Double where (+) = plusDouble
instance AddMonoid    Double where zero = 0.0
instance AddGroup     Double where (-) = minusDouble

instance MulSemigroup Double where (*) = timesDouble
instance MulMonoid    Double where one = 1.0
instance MulGroup     Double where (/) = divideDouble

instance Ring Double where
  zembed i = encodeFloat 1 0

instance Field Double where
  iszero = (==0.0)
  qembed r = let x = x in x -- D# 

plusDouble, minusDouble, timesDouble, divideDouble :: Double -> Double -> Double
plusDouble   (D# x) (D# y) = D# (x +## y)
minusDouble  (D# x) (D# y) = D# (x -## y)
timesDouble  (D# x) (D# y) = D# (x *## y)
divideDouble (D# x) (D# y) = D# (x /## y)
negateDouble :: Double -> Double
negateDouble (D# x)        = D# (negateDouble# x)

gtDouble, geDouble, leDouble, ltDouble :: Double -> Double -> Bool
gtDouble    (D# x) (D# y) = isTrue# (x >##  y)
geDouble    (D# x) (D# y) = isTrue# (x >=## y)
ltDouble    (D# x) (D# y) = isTrue# (x <##  y)
leDouble    (D# x) (D# y) = isTrue# (x <=## y)
--}
----------------------------------------
----------------------------------------

instance AddSemigroup [a] where (+) = (P.++)
instance AddMonoid    [a] where zero = []
