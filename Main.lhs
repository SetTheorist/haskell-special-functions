\documentclass{article}
\usepackage[top=1in,bottom=1in,left=1in,right=1in]{geometry}
\usepackage{amsmath}
\usepackage{amssymb}

%\usepackage{verbatim}
%\newenvironment{code}{\footnotesize\verbatim}{\endverbatim\normalsize}

%\usepackage{listings}
%\lstnewenvironment{code}{\lstset{language=Haskell,basicstyle=\small}}{}
%basicstyle=\small\ttfamily,
\usepackage{listings}
\lstloadlanguages{Haskell}
\lstnewenvironment{code}
    {\lstset{}%
      \csname lst@SetFirstLabel\endcsname}
    {\csname lst@SaveFirstLabel\endcsname}
    \lstset{
      language=Haskell,
      basicstyle=\small,
      flexiblecolumns=false,
      basewidth={0.5em,0.45em},
      literate={+}{{$+$}}1 {/}{{$/$}}1 {*}{{$*$}}1 {=}{{$=$}}1
               {>}{{$>$}}1 {<}{{$<$}}1 {\\}{{$\lambda$}}1
               {\\\\}{{\char`\\\char`\\}}1
               {->}{{$\rightarrow$}}2 {>=}{{$\geq$}}2 {<-}{{$\leftarrow$}}2
               {<=}{{$\leq$}}2 {=>}{{$\Rightarrow$}}2 
               {\ .}{{$\circ$}}2
               {>>}{{>>}}2 {>>=}{{>>=}}2
               %{\ .\ }{{$\circ$}}2
               %{|}{{$\mid$}}1               
               % my modifications:
               {\ .\ }{{$\circ$}}3
               {\$}{{\scriptsize\$}}1
               {!}{{\scriptsize!}}1
               {&&}{{$\wedge$}}2               
               {||}{{$\vee$}}2               
               %{'}{{${}^\prime$}}1 %% TODO: this is ugly, find better
    }

\DeclareMathOperator{\agm}{agm}
\DeclareMathOperator{\Ai}{Ai}
\DeclareMathOperator{\Bi}{Bi}
\DeclareMathOperator{\erf}{erf}
\DeclareMathOperator{\erfc}{erfc}
\DeclareMathOperator{\Li}{Li}

\title{Computation of Special Functions\\(Haskell)}
\author{Apollo Hogan}
\begin{document}

\maketitle
\tableofcontents

\section{Introduction}

Special functions.

%\section{Util}
\input{Util.lhs}

%\section{Fibo}
\input{Fibo.lhs}

%\section{Numbers}
\input{Numbers.lhs}

%\section{Exponential, logarithm}
\input{Exp.lhs}

%\section{Gamma}
\input{Gamma.lhs}

%\section{Erf}
\input{Erf.lhs}

%\section{AGM}
\input{AGM.lhs}

%\section{Airy}
\input{Airy.lhs}

%\section{Riemann Zeta}
\input{Zeta.lhs}

%\section{Spence}
\input{Spence.lhs}

%\section{Lommel}
\input{Lommel.lhs}

\end{document}
