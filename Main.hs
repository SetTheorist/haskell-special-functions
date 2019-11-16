import Control.Monad

import SpecialFunctions

main :: IO ()
main = do
    forM_ [-10,-9.99..5] $ \x -> do
        let ax = sf_airy_ai x
        let bx = sf_airy_bi x
        putStr $ (show x)++" "++(show ax)++" "++(show bx)++"\n"

{--
main :: IO ()
main = do
    forM_ [-20,-19.99..20] $ \x -> do
        let ex = sf_erf x
        let cx = sf_erfc x
        putStr $ (show x)++" "++(show ex)++" "++(show cx)++"\n"
--}
