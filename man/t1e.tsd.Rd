\encoding{utf-8}
\name{t1e.tsd}
\alias{t1e.tsd}
\title{Significance level adjustment for bioequivalence studies using adaptive two-stage 2x2 crossover designs}
\description{
This function calculates, by iterative search, the adjusted significance levels to be used in each stage of adaptive two-stage 2\enc{×}{x}2 crossover designs, to ensure an overall \ifelse{html}{\out{type&nbsp;I}}{type I} error below a specified significance level.
}
\usage{
t1e.tsd(n1, CV, GMR = 0.95, Nmax = 150, min.n2 = n1/2, type = 1,
        alpha = 0.05, alpha1 = 0.0294, alpha2 = alpha1, targetpower = 0.8,
        setseed = TRUE, details = TRUE, print = TRUE, ...)
}
\arguments{
\item{n1}{Sample size of \ifelse{html}{\out{stage&nbsp;1}}{stage 1}.}

\item{CV}{Within subject coeficient of variation (use \emph{e.g.}, 0.3 for 30\%).}

\item{GMR}{Expected geometric mean ratio to be used in decision scheme (power calculations in
    \ifelse{html}{\out{stage&nbsp;1}}{stage 1} and sample size estimation for
    \ifelse{html}{\out{stage&nbsp;2}}{stage 2}).\cr
    By default \verb{0.95}.}

\item{Nmax}{Overall maximum number of subjects (sum of sample sizes in both stages, \emph{i.e.}, \verb{N = n1 + N2}).\cr
By default \verb{150}, see the details section.}

\item{min.n2}{Minimum number of subjects at second stage. By default \verb{n1/2}. Set \verb{min.n2 = 0} to cancel any limitation on the sample size at the second stage, \verb{N2}. See the details section for more information.}

\item{type}{Type 1 or 2 methodology. See the details section. By default \verb{1}.}

\item{alpha}{Target overall significance level (both stages). By default \verb{0.05}.}

\item{alpha1}{Initial significance level at \ifelse{html}{\out{stage&nbsp;1}}{stage 1}. By default \verb{0.0294}.}

\item{alpha2}{Initial significance level at \ifelse{html}{\out{stage&nbsp;2}}{stage 2}. By default \verb{0.0294}.}

\item{targetpower}{Power threshold in the power monitoring steps and power to achieve in the sample size estimation step.\cr
By default \verb{0.8}.}

\item{setseed}{Simulations are dependent on the starting point of the (pseudo) random number generator.
To avoid differences in power for different runs a setseed(1234567) is issued if set to \verb{TRUE}, the default.\cr
Set this argument to \verb{FALSE} for a random seed.}

\item{details}{If set to \verb{TRUE} (default) shows intermediate results in the console. Set this argument to \verb{FALSE} to supress intermediate results.}

\item{print}{If set to \verb{TRUE} (default) shows final results in the console.\cr
Set this argument to \verb{FALSE} to return a list of final results.}

\item{...}{Optional additional arguments. See package Power2Stage, function \link[Power2Stage:power.tsd]{power.tsd}.}
}
\value{
Adjusted significance levels, \verb{alpha1} and \verb{alpha2}.
}
\details{
The type 1 method uses the same adjusted \verb{alpha1} and \verb{alpha2} in \ifelse{html}{\out{stages&nbsp;1}}{stages 1} and 2, respectively.\cr
The type 2 method uses an unadjusted \verb{alpha} if interim power is at least the targetpower, or adjusted \verb{alpha1} and \verb{alpha2} in \ifelse{html}{\out{stages&nbsp;1}}{stages 1} and 2 otherwise. The terminology follows \enc{Schütz}{Schuetz} (2015).

By default, \verb{Nmax = 150} and \verb{min.n2 = n1/2}. Its present implementation is based on controlling the total sample size: After computing the required sample size in the second step, say \verb{N2} (to ensure the required power \verb{targetpower}), the number of additional observations is computed as \verb{N2 = max(min.n2, N2)}, but if \verb{n1 + n2 > Nmax}, the study is terminated. These criteria and default values are based on Molins \emph{et al.} (2017).

\verb{n1} and \verb{CV} can be scalars or vectors. As vectors, the result might loose some power for a particular combination.\cr
Global significance level is ensured by returning the adjusted \verb{alpha1} and \verb{alpha2} of the worst case scenario for all combinations of \verb{n1} and \verb{CV}.

By default \verb{alpha = 0.05}, \verb{alpha1 = 0.0294} and \verb{alpha2 = alpha1}. However, \verb{alpha1} and \verb{alpha2} can be different.\cr
In this case, set \verb{alpha1} and \verb{alpha2}, and \verb{alpha1} will be fixed whilst \verb{alpha2} adjusted.

This function uses the package 'Power2Stage' developed by Labes \emph{et al.} and is based on the Potvin \emph{et al.} (2008) methods.
}
\examples{
\dontrun{
t1e.tsd(n1 = 24, CV = c(0.3, 0.4, 0.5, 0.6), GMR = 0.95, type = 1)
# should give adjusted alpha 0.0308 at both stages and
# a maximum empirical type I error 0.049807
# Note: runtime ~11 minutes on a Xeon E3-1245 Quadcore 3.4 GHz
t1e.tsd(n1 = 12, CV = 0.2, GMR = 0.9, alpha1 = 0.0280, alpha2 = 0.0280, type = 2)
# should give adjusted alpha 0.0269 at both stages and
# a maximum empirical type I error 0.04997
# Note: runtime ~45 minutes
# List of results
x <- t1e.tsd(n1 = 12, CV = 0.10, type = 2, print = FALSE)
print(as.data.frame(x), row.names = FALSE)}
}
\references{
Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs.}\cr
Pharm Stat. 2008; 7(4):245--62. \doi{10.1002/pst.294}

Montague TH, Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ. \emph{Additional results for \sQuote{Sequential design approaches for bioequivalence studies with crossover designs}.}\cr
Pharm Stat. 2011; 11(1):8--13. \doi{10.1002/pst.483}

Fuglsang A. \emph{Controlling type I errors for two-stage bioequivalence study designs.}\cr
Clin Res Reg Aff. 2011; 28(4):100--5. \doi{10.3109/10601333.2011.631547}

\enc{Schütz}{Schuetz} H. \emph{Two-stage designs in bioequivalence trials.}\cr
Eur J Clin Pharmacol. 2015; 71(3):271--81. \doi{10.1007/s00228-015-1806-2}

Molins E, Cobo E, \enc{Ocaña}{Ocana} J. \emph{Two-stage designs versus European scaled average designs in bioequivalence studies for highly variable drugs: Which to choose?}\cr
Stat Med. 2017; 36(30):4777--88. \doi{10.1002/sim.7452}
}
