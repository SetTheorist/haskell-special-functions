\chapter{Solving}

\section{Preamble}
\begin{code}
{-# Language BangPatterns #-}
module Solve where
import Trig
import Util
import System.IO.Unsafe
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Linear}
Solve the equation $a + bx = 0$
\begin{code}
sf_solve_linear :: (Value v) => v -> v -> v
sf_solve_linear !a !b = -a/b
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Quadratic}
Solve the equation $a + bx + cx^2 = 0$
\begin{code}
sf_solve_quadratic :: (Value v) => v -> v -> v -> [v]
sf_solve_quadratic !a !b !c
  | c==0 =
      let !r1 = sf_solve_linear a b
      in [r1]
  | otherwise =
  -- TODO: make this robust!
      let !r1 = (-b + (sf_sqrt $ b^2 - 4*a*c)) / (2*c)
          !r2 = (-b - (sf_sqrt $ b^2 - 4*a*c)) / (2*c)
      in [r1, r2]
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Cubic}
Solve the equation $dx^3 + cx^2 + bx + a = 0$
currently uses the trigonometric approach.
TODO: make robust!
  fails, for example, on $x^3 + 3x^2 + 3x + 1 = (x+1)^3$
  (get division-by-zero and NaN result)
  (should eliminate repeated roots first!)
\begin{code}
sf_solve_cubic :: (Value v) => v -> v -> v -> v -> [v]
sf_solve_cubic !a !b !c !d
  | d==0 = sf_solve_quadratic a b c
  | otherwise =
      -- get equivalent "depressed" cubic:
      --   t^3 + pt + q = 0
      let !p = (3*d*b - (c^2)) / (3*(d^2))
          !q = (2*(c^3) - 9*d*c*b + 27*(d^2)*a) / (27*(d^3))
          -- trigonometric approach
          !t = sf_sqrt $ -4*p/3
          !alpha = sf_acos $ -4*q/(t^3)
          !r1' = t*(sf_cos $ alpha/3)
          !r2' = t*(sf_cos $ alpha/3 + 2*pi/3)
          !r3' = t*(sf_cos $ alpha/3 + 4*pi/3)
          !r1 = r1' - c/(3*d)
          !r2 = r2' - c/(3*d)
          !r3 = r3' - c/(3*d)
      in [r1,r2,r3]
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Quartic}
Solve the equation $ex^4 + dx^3 + cx^2 + bx + a = 0$
approach via factoring into quadratics
TODO: make robust!
blows up for some cases, e.g. $(x+1)^4 = x^4+4x^3+6x^2+4x+1$
Approach from wikipedia "Quartic function" article
\begin{code}
sf_solve_quartic :: (Value v) => v -> v -> v -> v -> v -> [v]
sf_solve_quartic !a !b !c !d !e
  | e==0 = sf_solve_cubic a b c d
  | otherwise =
    let --normalize
        !a' = a/e
        !b' = b/e
        !c' = c/e
        !d' = d/e
        !e' = 1
        -- transform to "depressed" quartic
        --  x^4 + C x^2 + D x + E
        !cc = c' - 3*d'^2/8
        !dd = d'^3/8 - d'*c'/2 + b'
        !ee = a' - 3*d'^4/256 + d'^2*c'/16 - b'*d'/4

        --  solve resolvent cubic
        --  P^3 + 2cP^2 (c^2-4e)P - d^2 = 0
        --  (P=p^2)
        [p1,p2,p3] = sf_solve_cubic (-dd^2) (cc^2-4*ee) (2*cc) (1)
        -- err1 = P1^3 + (2*C)*P1^2 + (C^2-4*E)*P1 - D^2
        -- err2 = P2^3 + (2*C)*P2^2 + (C^2-4*E)*P2 - D^2
        -- err3 = P3^3 + (2*C)*P3^2 + (C^2-4*E)*P3 - D^2

        !p = sf_sqrt p3
        !r = -p
        !s = (cc + p^2 + dd/p)/2
        !q = (cc + p^2 - dd/p)/2

        -- thus x^4+bx^3+cx^2+dx+e = (x^2+px+q)(x^2+rx+s)
        [r1,r2] = sf_solve_quadratic q p 1
        [r3,r4] = sf_solve_quadratic s r 1

        !r1' = r1 - d/4
        !r2' = r2 - d/4
        !r3' = r3 - d/4
        !r4' = r4 - d/4
    in [r1',r2',r3',r4']
\end{code}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Polynomial}
Compute the root(s) of the polynomial $a_1 + a_2*x + ... + a_n*x^{n-1}$.
returns up to nr roots.
\begin{code}
sf_solve_poly :: (Value v) => Int -> [v] -> [v]
sf_solve_poly nr as = undefined
{--
## @deftypefn {Function File} {@var{res} =} sf_solve_poly (@var{a}, [@var{nr}])
## Compute the root(s) of the polynomial $a_1 + a_2*x + ... + a_n*x^(n-1)$.
## returns up to nr roots.
## @end deftypefn
function res = sf_solve_poly(a, nr)
  if (nargin < 1) print_usage; endif
  #if (nargin < 2) nr = 1; endif
  m = length(a)-1; # degree

  #nr = min(nr, m);
  #res = ones(nr);
  #res(1) = laguerre(a, 1);

  res = ones(m,1);
  ad = a;
  for jj = (m-1):(-1):0
    x = 0.0;
    ad_v = ad(1:(jj+2));
    x = laguerre(ad_v, x);
    if (abs(imag(x)) <= 1e-16*abs(real(x)))
      x = real(x);
    endif
    res(jj+1) = x;
    # forward deflation:
    b = ad(jj+2);
    for jjj = jj:(-1):0
      c = ad(jjj+1);
      ad(jjj+1) = b;
      b = x*b + c;
    endfor
  endfor

  # polish the roots
  da = (1:m) .* a(2:end);
  for jj = 1:m
    res(jj) = polish(a, da, res(jj));
  endfor
  res = sort(res);
endfunction

function res = polish(a, da, x)
  # try one Newton step, else Laguerre with full coeffs
  fx = sf_polynomial_value(a, x);
  dfx = sf_polynomial_value(da, x);
  new_x = x - fx/dfx;
  f_newx = sf_polynomial_value(a, new_x);
  if (abs(f_newx) < abs(fx))
    res = new_x;
  else
    res = laguerre(a, x);
  endif
endfunction

# Laguerre's method, adapted from NR
function res = laguerre(a, x)
  persistent MR=8;
  persistent MT=12;
  persistent MAXIT=MR*MT;
  persistent frac = [0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0];
  persistent eps = 1e-17;
  m = length(a);
  for it = 1:MAXIT
    b = a(m);
    err = abs(b);
    d = f = 0;
    abx = abs(x);
    for jj = m-1 : (-1) : 1
      f = x*f + d;
      d = x*d + b;
      b = x*b + a(jj);
      err = abs(b) + abx*err;
    endfor
    err *= eps;
    if (abs(b) <= err)
      #on a root
      break;
    endif
    g = d/b;
    g2 = g*g;
    h = g2 - 2*f/b;
    sq = sqrt((m-1)*m*h - g2);
    gp = g + sq;
    gm = g - sq;
    abp = abs(gp);
    abm = abs(gm);
    if (abp < abm) gp = gm; endif
    if (max(abp,abm) > 0)
      dx = m/gp;
    else
      dx = (1+abx)*exp(I*it);
    endif
    x1 = x - dx;
    if (x == x1)
      #converged
      break;
    endif
    if (rem(it, MT)!=0) x = x1;
    else x -= frac(it/MT) * dx;
    endif
  endfor
  res = x;
endfunction
--}
\end{code}
