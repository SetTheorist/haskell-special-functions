\section{AGM}

\subsection{Preamble}

\begin{titled-frame}{\color{blue}\tt module AGM}%
\begin{code}
module AGM (sf_agm, sf_agm') where
import Util
\end{code}
\end{titled-frame}

\subsection{AGM}
Gauss' arithmetic-geometric mean or AGM of two numbers is defined as
the limit $\agm(\alpha,\beta)=\lim_n\alpha_n=\lim_n \beta_n$ where
we define
\begin{eqnarray*}
  \alpha_{n+1} &=& \frac{\alpha_n+\beta_n}{2} \\
  \beta_{n+1} &=& \sqrt{\alpha_n\cdot\beta_n}
\end{eqnarray*}
(Note that we need real values to be positive for this to make sense.)

\subsubsection{\tt sf\_agm alpha beta}
Here we compute the AGM via the definition and return the full
arrays of intermediate values $([\alpha_n],[\beta_n],[\gamma_n])$,
where $\gamma_n=\frac{\alpha_n-\beta_n}{2}$.
(The iteration converges quadratically so this is an efficient approach.)
\begin{titled-frame}{$\text{\tt sf\_agm alpha beta}=\agm(\alpha,\beta)$\marginnote{\tt sf\_agm}}%
\begin{code}
sf_agm :: (Value v) => v -> v -> ([v],[v],[v])
sf_agm alpha beta = agm [alpha] [beta] [alpha-beta]
  where agm as@(a:_) bs@(b:_) cs@(c:_) = 
          if c==0 then (as,bs,cs)
          else let a' = (a+b)/2
                   b' = sf_sqrt (a*b)
                   c' = (a-b)/2
               in if c'==c then (as,bs,cs)
                  else agm (a':as) (b':bs) (c':cs)
\end{code}%
\end{titled-frame}
  
\subsubsection{\tt sf\_agm' alpha beta}
Here we return simply the value $\verb|sf_agm' a b|=\agm(a,b)$.
\begin{titled-frame}{$\text{\tt sf\_agm' z}=\agm z$}%
\begin{code}
sf_agm' :: (Value v) => v -> v -> v
sf_agm' alpha beta = agm alpha beta ((alpha-beta)/2)
  --let (as,_,_) = sf_agm alpha beta in head as
  where agm a b 0 = a
        agm a b c =
          let a' = (a+b)/2
              b' = sf_sqrt (a*b)
              c' = (a-b)/2
          in agm a' b' c'
\end{code}
\end{titled-frame}

\begin{code}
sf_agm_c0 :: (Value v) => v -> v -> v -> ([v],[v],[v])
sf_agm_c0 alpha beta c0 = undefined
\end{code}
