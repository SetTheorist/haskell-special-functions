\chapter{Gamma}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Preamble}
A basic preamble.
\begin{titled-frame}{\color{blue}\tt module Gamma}
\begin{code}
module Gamma (
    euler_gamma,
    factorial,
    sf_beta,
    sf_gamma,
    sf_invgamma,
    sf_lngamma,
    sf_digamma,
    )
where
import Exp
import Numbers(factorial,sf_bernoulli_b)
import Trig
import Util
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Misc}

\subsection{\tt euler\_gamma}
A constant for Euler's gamma:
\[ \gamma = \lim_{n\to\infty}\left(\sum_{k=1}^n\frac1n - \ln n\right) \marginnote{$\gamma$}\]
\begin{code}
euler_gamma :: (Floating a) => a
euler_gamma = 0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749
\end{code}

\subsection{\tt sf\_beta a b}
The Beta integral
\[ B(a,b) = \int_0^1 t^{a-1}(1-t)^{b-1}\,dt = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} \marginnote{$B(a,b)$}\]
implemented in terms of log-gamma
\[ \verb|sf_beta a b| = e^{\ln\Gamma(a) + \ln\Gamma(b) - \ln\Gamma(a+b)} \]
\begin{code}
sf_beta :: (Value v) => v -> v -> v
sf_beta a b = sf_exp $ (sf_lngamma a) + (sf_lngamma b) - (sf_lngamma$a+b)
\end{code}



Compute the Pochhammer symbol $(z)_\alpha$
\begin{code}
{--
sf_pochhammer z al
  | (sf_is_nonnegint al) && al<100 = product $ map (z+) (0:(floor$al-1))
  | otherwise = sf_exp $ (sf_lngamma $ z + al) - (sf_lngamma z)
  --  kludge cleanup
  -- f (sf_is_int(z) && sf_is_int(al)) res = round(res); endif
  -- f (imag(z)==0 && imag(al)==0) res = real(res); endif
--}
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Gamma}

The gamma function
\[ \Gamma(z) = \int_0^\infty e^{-t}t^{z} \,\frac{dt}{t} \marginnote{$\Gamma(z)$}\]

\subsection{\tt sf\_gamma z}
The gamma function implemented using the identity $\Gamma(z) = \frac{1}{z}\Gamma(z+1)$
to increase the real part of the argument to be $>15$ and then
using an asymptotic expansion for log-gamma, \verb|lngamma__asymp|, to evaluate.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_gamma x} = \Gamma(x)$}
\begin{code}
sf_gamma :: (Value v) => v -> v
sf_gamma x =
  redup x 1 $ \ x' t -> t * (sf_exp (lngamma__asymp x'))
  where redup x t k
          | (re x)>15 = k x t
          | otherwise = redup (x+1) (t/x) k
\end{code}
\end{titled-frame}

\subsubsection{\tt lngamma\_\_asymp z}
The asymptotic expansion for log-gamma
\[ \ln\Gamma(z) \sim (z-\frac12)\ln z - z + \frac12\ln(2\pi) + \sum_{k=1}^{\infty}\frac{B_{2k}}{2k(2k-1)z^{2k-1}} \]
where $B_n$ is the $n$'th Bernoulli number.
\begin{titled-frame}{\tt lngamma\_\_asymp z}
\begin{code}
lngamma__asymp :: (Value v) => v -> v
lngamma__asymp z = (z - 1/2)*(sf_log z) - z + (1/2)*sf_log(2*pi) + (ksum terms)
  where terms = [b2k/(2*k*(2*k-1)*z^(2*k'-1)) | k'<-[1..10], let k=(#)k', let b2k=fromRational$sf_bernoulli_b!!(2*k')]
\end{code}
\end{titled-frame}

\subsection{\tt sf\_invgamma z}
The inverse gamma function, $\verb|sf_invgamma z|=\frac{1}{\Gamma(z)}$.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_invgamma x} = 1/\Gamma(x)$}
\begin{code}
sf_invgamma :: (Value v) => v -> v
sf_invgamma x =
  let (x',t) = redup x 1
      lngx = lngamma__asymp x'
  in t * (sf_exp$ -lngx)
  where redup x t
          | (re x)>15 = (x,t)
          | otherwise = redup (x+1) (t*x)
\end{code}
\end{titled-frame}

\subsection{\tt sf\_lngamma z}
The log-gamma function, $\verb|sf_lngamma z|=\ln\Gamma(z)$.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_lngamma x} = \ln\Gamma(x)$}
\begin{code}
sf_lngamma :: (Value v) => v -> v
sf_lngamma x =
  let (x',t) = redup x 0
      lngx = lngamma__asymp x'
  in t + lngx
  where redup x t
          | (re x)>15 = (x,t)
          | otherwise = redup (x+1) (t-sf_log x)

\end{code}
\end{titled-frame}

\subsubsection{Spouge's approximation to the gamma function}
In tests, this gave disappointing results.
\begin{code}
-- Spouge's approximation (a=17?)
spouge_approx :: (Value v) => Int -> v -> v
spouge_approx a z' =
  let z = z' - 1
      a' = (#)a
      res = (z+a')**(z+(1/2)) * sf_exp (-(z+a'))
      sm = fromDouble$sf_sqrt(2*pi)
      terms = [(spouge_c k a') / (z+k') | k<-[1..(a-1)], let k' = (#)k]
      smm = sm + ksum terms
  in res*smm
  where
    spouge_c k a = ((if k`mod`2==0 then -1 else 1) / ((#) $ factorial (k-1)))
                    * (a-((#)k))**(((#)k)-1/2) * sf_exp(a-((#)k))
\end{code}

\begin{code}
spouge :: (Value v) => Int -> v -> v
spouge a' z' =
  let z = z' - 1
      a = fromDouble$(#)a'
      -- I don't quite understand why I can't do this:
      --q = fromReal $ (sf_sqrt(2*pi) :: (RealKind v))
      q = sf_sqrt(2*pi)
  in (z+a)**(z+1/2)*(sf_exp(-z-a))*(q + ksum (map (\k->(c a k)/(z+(#)k)) [1..(a'-1)]))
  where
    c :: (Value v) => v -> Int -> v
    c a k = let k' = (#)k
                sgn = if even k then -1 else 1
            in sgn*(a-k')**(k'-1/2)*(sf_exp(a-k')) / ((#)$factorial(k-1))
\end{code}  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Digamma}

The digamma function
\[ \psi(z) = \frac{d}{dz} \ln\Gamma(z) = \frac{\Gamma'(z)}{\Gamma(z)} \marginnote{$\psi(z)$}\]

\subsection{\tt sf\_digamma z}
We implement with a series expansion for $|z|<=10$ and otherwise with an asymptotic expansion.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_digamma z} = \psi(z)$}
\begin{code}
sf_digamma :: (Value v) => v -> v
--sf_digamma n | is_nonposint n = Inf
sf_digamma z | (rabs z)>10 = digamma__asympt z
             | otherwise   = digamma__series z
\end{code}
\end{titled-frame}

\subsubsection{\tt digamma\_\_series z}
The series expansion is the following
\[ \psi(z) = -\gamma-\frac1z + \sum_{k=1}^\infty \frac{z}{k(k+z)} \]
but with Euler-Maclaurin correction terms:
\[ \psi(z) = -\gamma-\frac1z + \sum_{k=1}^n\frac{z}{k(k+z)}
    + (\ln\frac{k+z}{k} - \frac{z}{2k(k=z)} + \sum_{j=1}^{p}B_{2j}(k^{-2j}-(k+z)^{-2j}) \]
\begin{titled-frame}{\color{blue}\tt digamma\_\_series z}
\begin{code}
digamma__series :: (Value v) => v -> v
digamma__series z =
  let res = -euler_gamma - (1/z)
      terms = map (\k->z/((#)k*(z+(#)k))) [1..]
      corrs = map (correction.(#)) [1..]
  in summer res res terms corrs
  where
    summer :: (Value v) => v -> v -> [v] -> [v] -> v
    summer res sum (t:terms) (c:corrs) =
      let sum' = sum + t
          res' = sum' + c
      in if res==res' then res
         else summer res' sum' terms corrs
    bn1 = fromRational$sf_bernoulli_b!!2
    bn2 = fromRational$sf_bernoulli_b!!4
    bn3 = fromRational$sf_bernoulli_b!!6
    bn4 = fromRational$sf_bernoulli_b!!8
    correction k =
      (sf_log$(k+z)/k) - z/2/(k*(k+z))
        + bn1*(k^^(-2) - (k+z)^^(-2))
        + bn2*(k^^(-4) - (k+z)^^(-4))
        + bn3*(k^^(-6) - (k+z)^^(-6))
        + bn4*(k^^(-8) - (k+z)^^(-8))
\end{code}
\end{titled-frame}

\subsubsection{\tt digamma\_\_asympt z}
The asymptotic expansion (valid for $|arg z|<\pi$) is the following
\[ \psi(z) \sim \ln z - \frac{1}{2z} + \sum_{k=1}^\infty\frac{B_{2k}}{2kz^{2k}} \]
If $\Re z<\frac12$ then we use the reflection identity to ensure $\Re z\geq\frac12$:
\[ \psi(z) - \psi(1-z) = \frac{-\pi}{\tan(\pi z)} \]
\begin{titled-frame}{\color{blue}\tt digamma\_\_asympt z}
\begin{code}
digamma__asympt :: (Value v) => v -> v
digamma__asympt z
  | (re z)<0.5 = compute (1 - z) $ -pi/(sf_tan(pi*z)) + (sf_log(1-z)) - 1/(2*(1-z))
  | otherwise  = compute z $ (sf_log z) - 1/(2*z)
    where
      compute z res =
        let z_2 = z^^(-2)
            zs = iterate (*z_2) z_2
            terms = zipWith (\n z2n -> z2n*(fromRational$sf_bernoulli_b!!(2*n+2))/(#)(2*n+2)) [0..] zs
        in sumit res res terms
      sumit res ot (t:terms) =
        let res' = res - t
        in if res==res' || (rabs t)>(rabs ot)
           then res
           else sumit res' t terms
\end{code}
\end{titled-frame}

