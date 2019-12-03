%\documentclass{article}
\documentclass{report}
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\F}[1]{{\mathtt{#1}}}
\DeclareMathOperator{\erf}{erf}
\DeclareMathOperator{\erfc}{erfc}
\DeclareMathOperator{\acos}{acos}
\DeclareMathOperator{\asin}{asin}
\DeclareMathOperator{\atan}{atan}
\DeclareMathOperator{\acosh}{acosh}
\DeclareMathOperator{\asinh}{asinh}
\DeclareMathOperator{\atanh}{atanh}
\DeclareMathOperator{\ph}{ph}
\let\O=\operatorname
\newcommand{\Ob}[1]{\operatorname{\mathbf{#1}}}
\newcommand{\RR}{{\mathbb{R}}}
\newcommand{\CC}{{\mathbb{C}}}
\newcommand{\ZZ}{{\mathbb{Z}}}
\newcommand{\QQ}{{\mathbb{Q}}}
\newcommand{\NN}{{\mathbb{N}}}
\newcommand{\NNo}{{\mathbb{N}_0}}
\newcommand{\ii}{{\hat{\imath}}}
\newcommand{\Hyper}[5]{{\mathop{{}_{#1}\!\O{F\!}_{#2}\left({\genfrac{}{}{0pt}{}{#3}{#4}}\middle|{#5}\right)}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\title{Computation of Special Functions\\(Haskell)}
\author{Apollo Hogan}
\begin{document}

\maketitle
\tableofcontents

\chapter{Introduction}

Special functions.

%\chapter{Util}
\input{Util.lhs}

%\chapter{Fibo}
\input{Fibo.lhs}

%\chapter{Numbers}
\input{Numbers.lhs}

%\chapter{Exponential, logarithm}
\input{Exp.lhs}

%\chapter{Gamma}
\input{Gamma.lhs}

%\chapter{IncompleteGamma}
\input{IncompleteGamma.lhs}

%\chapter{Erf}
\input{Erf.lhs}

%\chapter{Bessel}
\input{Bessel.lhs}

%\chapter{ExpInt}
\input{ExpInt.lhs}

%\chapter{AGM}
\input{AGM.lhs}

%\chapter{Airy}
\input{Airy.lhs}

%\chapter{Riemann Zeta}
\input{Zeta.lhs}

%\chapter{Elliptic}
\input{Elliptic.lhs}

%\chapter{Debye}
\input{Debye.lhs}

%\chapter{Spence}
\input{Spence.lhs}

%\chapter{Lommel}
\input{Lommel.lhs}

\end{document}
