module OrthoPoly where

class OrthogonalPolynomial a where
  in_support :: a -> Double -> Bool
  weight :: a -> Double -> Double

  coeffs :: a -> Int -> [Double]
  scale :: a -> Int -> Double
  value :: a -> Int -> Double -> Double
  weights :: a -> Int -> [Double]
  zeros :: a -> Int -> [Double]

  value_arr :: a -> Int -> Double -> [Double]

----------------------------------------

data LegendrePolynomial = LegendrePolynomial

instance OrthogonalPolynomial LegendrePolynomial where
  in_support :: a -> Double -> Bool
  weight :: a -> Double -> Double

  coeffs :: a -> Int -> [Double]
  scale :: a -> Int -> Double
  value :: a -> Int -> Double -> Double
  weights :: a -> Int -> [Double]
  zeros :: a -> Int -> [Double]

  value_arr :: a -> Int -> Double -> [Double]

----------------------------------------

