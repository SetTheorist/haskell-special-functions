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
    }

\begin{document}

\section{Introduction}

Special functions.

%\section{Util}
\input{Util.lhs}

%\section{Fibo}
\input{Fibo.lhs}

%\section{Exponential, logarithm}
\input{Exp.lhs}

\end{document}
