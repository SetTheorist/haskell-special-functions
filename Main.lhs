\documentclass{article}
\usepackage[top=1in,bottom=1in,left=1in,right=1in]{geometry}
\usepackage{amsmath}
\usepackage{amssymb}
%\usepackage{mathbbol}

\usepackage{marginnote}

\usepackage{xcolor}
\usepackage{framed}
% to use titled-frame, must define TFFrameColor and TFTitleColor
\definecolor{TFFrameColor}{rgb}{0.9,0.9,0.9}
\definecolor{TFTitleColor}{rgb}{0.0,0.0,0.0}

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
               {>}{{$>$}}1 {<}{{$<$}}1 {\\}{{${\boldsymbol\lambda}$}}1
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
               {forall}{{$\forall$}}1%6
               %{'}{{${}^\prime$}}1 %% TODO: this is ugly, find better
    }

\DeclareMathOperator{\agm}{agm}
\DeclareMathOperator{\Ai}{Ai}
\DeclareMathOperator{\Bi}{Bi}
\DeclareMathOperator{\Ei}{Ei}
\DeclareMathOperator{\erf}{erf}
\DeclareMathOperator{\erfc}{erfc}
\DeclareMathOperator{\gud}{gud}
\DeclareMathOperator{\Li}{Li}
\DeclareMathOperator{\ph}{ph}

%\def\ii{{\mathbb{i}}}
\def\ii{{\hat\imath}}

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

%\section{IncompleteGamma}
\input{IncompleteGamma.lhs}

%\section{Erf}
\input{Erf.lhs}

%\section{Bessel}
\input{Bessel.lhs}

%\section{ExpInt}
\input{ExpInt.lhs}

%\section{AGM}
\input{AGM.lhs}

%\section{Airy}
\input{Airy.lhs}

%\section{Riemann Zeta}
\input{Zeta.lhs}

%\section{Elliptic}
\input{Elliptic.lhs}

%\section{Spence}
\input{Spence.lhs}

%\section{Lommel}
\input{Lommel.lhs}

\end{document}
